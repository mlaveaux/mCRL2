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

/// \brief Print the intermediate matches that succeeded.
constexpr bool PrintMatchSteps = false;

/// \brief Print the intermediate steps performed during construction.
constexpr bool PrintConstructionSteps = false;

using namespace mcrl2::data;
using namespace mcrl2::data::detail;

using namespace mcrl2::log;

/// \brief Compute the set of positions of t such that t[p] is a variable for all p in fringe.
void fringe_impl(const atermpp::aterm_appl& appl, position current, std::set<position>& fringe)
{
  if (is_variable(appl))
  {
    // The current position has a variable.
    fringe.insert(current);
  }
  else
  {
    // Extend the position to be one deeper into the subterm.
    current.emplace_back(0);
    for (const atermpp::aterm& argument : appl)
    {
      fringe_impl(static_cast<const atermpp::aterm_appl&>(argument), current, fringe);
      ++current.back();
    }
  }
}

/// \brief Compute the set of positions of t such that t[p] is a variable for all p in fringe(t).
std::set<position> fringe(const atermpp::aterm_appl& appl)
{
  std::set<position> result;
  fringe_impl(appl, position(), result);
  return result;
}

/// \brief Decides whether the left and right terms unify, assuming that vars(left) and vars(right) are disjoint.
bool unify(const atermpp::aterm_appl& left, const atermpp::aterm_appl& right)
{
  if (left == right)
  {
    // If both sides are equivalent, then they match under any substitution. If one contains
    // a variable then that variable is also not bound in the other term.
    return true;
  }
  else if (is_variable(left) || is_variable(right))
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

/// \brief Returns t[pos], i.e., the term at the given position using a index to keep track of the pos.
std::optional<data_expression> at_position_impl(const data_expression& t, const position& pos, std::size_t index)
{
  if (pos.empty())
  {
    // t[emptypos] = t
    return t;
  }
  else
  {
    std::size_t arg = pos[index];

    if (arg < t.size())
    {
      const data_expression& u = static_cast<const data_expression&>(static_cast<atermpp::aterm_appl>(t)[arg]);
      if (pos.size() == index + 1)
      {
        return u;
      }
      else
      {
        return at_position_impl(u, pos, index + 1);
      }
    }
  }

  return {};
}

/// \brief Returns t[pos], i.e., the term at the given position.
std::optional<data_expression> at_position(const data_expression& t, const position& pos)
{
  return at_position_impl(t, pos, 0);
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

/// \brief Replace the position variable at the given position by the expression c.
data_expression assign_at_position(const data_expression& term, const position& pos, const data_expression& c)
{
  mutable_indexed_substitution<variable, data_expression> sigma;
  sigma[position_variable(pos)] = c;
  return replace_variables(term, sigma);
}

/// \returns True iff there exists l in L : (exists pos' <= pos : head(l[pos']) in V
bool has_variable_higher_impl(const data_expression& t, const position& pos, std::size_t index)
{
  // These two conditions check pos' <= pos.
  if (pos.empty() || index >= pos.size())
  {
    return false;
  }
  else if (is_variable(static_cast<data_expression>(t[pos[index]])))
  {
    return true;
  }
  else
  {
    assert(index < t.size());
    return has_variable_higher_impl(static_cast<data_expression>(t[pos[index]]), pos, index + 1);
  }
}

/// \returns True iff (exists pos' <= pos : head(l[pos']) in V
bool has_variable_higher(const application& appl, const position& pos)
{
  return has_variable_higher_impl(appl, pos, 0);
}

/// \returns The number of arguments for the given data expressions.
inline std::size_t get_arity(const data_expression& t)
{
  if (!is_application(t) || is_not_equal(t))
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
  if (!is_application(t) || is_not_equal(t))
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
  mcrl2::utilities::stopwatch construction;
  enumerator_identifier_generator generator("@");

  // Preprocess the term rewrite system.
  for (const data_equation& old_equation : equations)
  {
    // Rename the variables in the equation
    auto [equation, partition] = rename_variables_unique(make_linear(old_equation, generator));

    // Add the index of the equation
    m_linear_equations.emplace_back(linear_data_equation(equation, partition));
  }

  // Determine the index of not_equal.
  m_not_equal_index = mcrl2::core::index_traits<mcrl2::data::function_symbol, function_symbol_key_type, 2>::index(static_cast<const function_symbol&>(not_equal()));

  // Construct the automaton.
  try
  {
    construct_apma(m_automaton.root(), position_variable(position()));

    mCRL2log(info) << "EAPMA (states: " << m_automaton.states() << ", transitions: " << m_automaton.transitions() << ") construction took " << construction.time() << " milliseconds.\n"
      << " there are " << m_positions.size() << " positions indexed.\n";

    // Ensure that the subterm indexing can store all possible terms in the required places.
    m_subterms.resize(m_positions.size());
  }
  catch (const std::exception& ex)
  {
    mCRL2log(info) << "Construction failed with " << ex.what() << ".\n";
  }
}

template<typename Substitution>
void AdaptiveMatcher<Substitution>::match(const data_expression& term)
{
  m_matching_sigma.clear();
  m_match_set = nullptr;
  m_match_index = 0;

  // Start with the root state.
  std::size_t s = m_automaton.root();

  // The empty position is the first index.
  m_subterms[0] = term;

  while (true)
  {
    // Store the labelling of the current state for convenience.
    const apma_state& state = m_automaton.label(s);

    // If L(s) = L' for some pattern set L'.
    if (state.match_set.size() > 0)
    {
      break;
    }

    // t[p] is given by the subterm table at the current position.
    const data_expression& subterm = m_subterms[state.position];
    if (subterm == data_expression())
    {
      // head(t[p]) is not defined so matching is finished.
      return;
    }

    if (PrintMatchSteps) { mCRL2log(info) << "Matching m_subterms[" << state.position << "] = " << subterm << "\n"; }

    // Update the (position) variable for the current subterm.
    if (state.variable.defined())
    {
      m_matching_sigma[state.variable] = subterm;
    }

    bool found_transition = false;
    if (is_application(subterm))
    {
      const auto& appl = static_cast<const application&>(subterm);

      // If delta(s, f) = (s', update)
      auto [s_prime, found] = m_automaton.transition(s, get_head_index(appl.head()));
      if (found)
      {
        if (PrintMatchSteps) { mCRL2log(info) << "Took transition from " << s << " to " << s_prime << " with label " << appl.head() << "\n"; }

        found_transition = true;
        s = s_prime;

        // Insert the subterms onto the position mapping.
        auto it = state.argument_positions.begin();
        for (const atermpp::aterm& argument : appl)
        {
          assert(it != state.argument_positions.end());
          assert(*it < m_subterms.size());
          m_subterms[*it] = static_cast<const data_expression&>(argument);
          ++it;
        }
      }
    }
    else if (is_function_symbol(subterm))
    {
      // If delta(s0, a) is defined for some term a followed by suffix t'.
      auto [s_prime, found] = m_automaton.transition(s, get_head_index(subterm));
      if (found)
      {
        if (PrintMatchSteps) { mCRL2log(info) << "Took transition " << s << " to " << s_prime << " with label " << static_cast<function_symbol>(subterm) << "\n"; }

        found_transition = true;
        s = s_prime;
      }
    }

    if (!found_transition)
    {
      // If delta(s, \neq)
      auto [s_prime, found] = m_automaton.transition(s, m_not_equal_index);
      if (found)
      {
        if (PrintMatchSteps)
        {
          mCRL2log(info) << "Took transition from " << s << " to " << s_prime << " with label not_equal.\n";
        }

        s = s_prime;
      }

      if (PrintMatchSteps) { mCRL2log(info) << "Matching failed.\n"; }
      return;
    }

    // Reset this position as we will not inspect it again.
    m_subterms[state.position] = data_expression();
  }

  if (PrintMatchSteps) { mCRL2log(info) << "Matching succeeded.\n"; }

  // for p in P do sigma := sigma[x_p -> t[p]]
  for (const auto& [var, pos] : m_automaton.label(s).variables)
  {
    assert(m_subterms[pos] != data_expression());
    m_matching_sigma[var] = m_subterms[pos];
  }

  m_match_set = &m_automaton.label(s).match_set;
}

template<typename Substitution>
const extended_data_equation* AdaptiveMatcher<Substitution>::next(Substitution& matching_sigma)
{
  if (m_match_set != nullptr)
  {
    while (m_match_index < m_match_set->size())
    {
      matching_sigma = m_matching_sigma;
      const linear_data_equation& result = (*m_match_set)[m_match_index];
      ++m_match_index;

      if (!is_consistent(result.partition(), matching_sigma))
      {
        // This rule matched, but its variables are not consistent w.r.t. the substitution.
        continue;
      }

      return &result;
    }
  }

  return nullptr;
}

// Private functions

template<typename Substitution>
void AdaptiveMatcher<Substitution>::construct_apma(std::size_t s, data_expression pref)
{
  if (PrintConstructionSteps) { mCRL2log(info) << "state = " << s << "("; }

  // L := {l in L | l unifies with pref}
  std::vector<std::reference_wrapper<const linear_data_equation>> L;
  for (const linear_data_equation& equation : m_linear_equations)
  {
    if (unify(equation.equation().lhs(), pref))
    {
      L.emplace_back(equation);
    }
  }

  // intersection := intersection_{l in L} fringe(l)
  std::set<position> intersection;
  bool first = true;
  for (const linear_data_equation& equation : L)
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

  // F := fringe(pref) \ intersection_{l in L} fringe(l)
  std::set<position> F;
  std::set<position> fringe_pref = fringe(pref);
  std::set_difference(fringe_pref.begin(), fringe_pref.end(), intersection.begin(), intersection.end(), std::inserter(F, F.begin()));

  // The labelling for the current state.
  apma_state& state = m_automaton.label(s);

  // if F = emptyset
  if (F.empty())
  {
    // M := M[L := L[s -> L]

    // Postprocessing: R := {r_i | l_i in L(s)}
    state.match_set.insert(state.match_set.begin(), L.begin(), L.end());

    // Postprocessing: P := union r_i in R : fringe(r_i) \ {L(s') in path(s)}
    // This is equivalent? to (union r_i in R : vars(r_i)) intersection vars(pref).
    std::set<variable> vars_in_rhs;
    std::for_each(L.begin(), L.end(),
      [&](const linear_data_equation& equation)
      {
        std::set<variable> vars = data::find_free_variables(equation.equation().rhs());
        vars_in_rhs.insert(vars.begin(), vars.end());
      });

    // Convert the resulting set to a vector for fast evaluation.
    if (PrintConstructionSteps) { mCRL2log(info) << "P = "; }
    std::set<position> F = fringe(pref);
    for (const position& pos : F)
    {
      if (PrintConstructionSteps) { mCRL2log(info) << pos << ", "; }
      std::optional<data_expression> var = at_position(pref, pos);
      assert(var);
      assert(is_variable(var.value()));

      state.variables.emplace_back(std::make_pair(static_cast<variable>(var.value()), m_positions.insert(pos).first));
    }

    // L' := L'[s -> (R, P)]
    if (PrintConstructionSteps) { mCRL2log(info) << ") \n"; }
  }
  else
  {
    // pos := select(F). For now the left-most position (left-to-right automaton).
    position pos = *std::min_element(F.begin(), F.end(), less_than);

    // M := M[L := L[s -> pos]]
    state.position = m_positions.insert(pos).first;
    if (PrintConstructionSteps) { mCRL2log(info) << "L = " << pos << ")\n"; }

    // for f in F s.t. exist l in L : head(l[pos]) = f
    std::set<std::pair<mcrl2::data::function_symbol, std::size_t>> symbols;

    for (const linear_data_equation& equation : L)
    {
      // t := l[pos]
      std::optional<data_expression> t = at_position(static_cast<const application&>(equation.equation().lhs()), pos);

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
      std::size_t s_prime = m_automaton.add_state();

      // M := M[delta := (delta := delta(s, f) -> s')
      m_automaton.add_transition(s, mcrl2::core::index_traits<mcrl2::data::function_symbol, function_symbol_key_type, 2>::index(symbol), s_prime);

      if (PrintConstructionSteps) { mCRL2log(info) << "added transition from " << s << " with label " << static_cast<atermpp::aterm>(symbol) << " to state " << s_prime << ".\n"; }

      // pref := f(omega^ar(f)), but we store in the correct position variables directly in the prefix.
      std::vector<data_expression> arguments;

      if (arity > 0)
      {
        // In this case there are multiple variables are introduced, their positions are derived from current position.i where i is their
        // position in the application f(x_0, ..., x_n).
        position var_position = pos;
        var_position.push_back(0);
        for (std::size_t index = 0; index < arity; ++index)
        {
          var_position.back() = index + 1;
          arguments.push_back(position_variable(var_position));
        }

        construct_apma(s_prime, assign_at_position(pref, pos, application(symbol, arguments.begin(), arguments.end())));
        max_arity = std::max(max_arity, arity);
      }
      else
      {
        // In this case the symbol is just a function symbol.
        construct_apma(s_prime, assign_at_position(pref, pos, symbol));
      }
    }

    // Here, we also ensure that the arguments of this symbol can be stored in the subterm table during evaluation.
    position var_position = pos;
    var_position.push_back(0);
    for (std::size_t index = 0; index < max_arity; ++index)
    {
      var_position.back() = index + 1;
      state.argument_positions.emplace_back(m_positions.insert(var_position).first);
    }

    // if exists l in L : (exists pos' <= pos : head(l[pos']) in V then
    if (std::find_if(L.begin(),
      L.end(),
      [&](const linear_data_equation& equation)
      {
        return has_variable_higher(static_cast<application>(equation.equation().lhs()), pos);
      }) != L.end())
    {
      // M := M[S := (S cup {s'})], where s' is a fresh state.
      std::size_t s_prime = m_automaton.add_state();

      // M := M[delta := (delta := delta(s, \neq) -> s')
      m_automaton.add_transition(s, mcrl2::core::index_traits<mcrl2::data::function_symbol, function_symbol_key_type, 2>::index(not_equal()), s_prime);

      if (PrintConstructionSteps) { mCRL2log(info) << "added transition from " << s << " with label " << static_cast<atermpp::aterm>(not_equal()) << " to state " << s_prime << ".\n"; }

      construct_apma(s_prime, assign_at_position(pref, pos, not_equal()));
    }
  }
}

// Explicit instantiations.

template class mcrl2::data::detail::AdaptiveMatcher<mutable_indexed_substitution<>>;
