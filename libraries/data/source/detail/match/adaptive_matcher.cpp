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
constexpr bool PrintMatchSteps = true;

/// \brief Print the intermediate steps performed during construction.
constexpr bool PrintConstructionSteps = true;

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
std::optional<data_expression> at_position_impl(const application& appl, const position& pos, std::size_t index)
{
  if (pos.empty())
  {
    // t[emptypos] = t
    return appl;
  }
  else
  {
    std::size_t arg = pos[index];

    if (arg < appl.size())
    {
      if (pos.size() == index + 1)
      {
        return appl[arg - 1];
      }
      else
      {
        return at_position_impl(static_cast<const application&>(appl[arg]), pos, index + 1);
      }
    }
  }

  return {};
}

/// \brief Returns t[pos], i.e., the term at the given position.
std::optional<data_expression> at_position(const application& appl, const position& pos)
{
  return at_position_impl(appl, pos, 0);
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
bool has_variable_higher_impl(const application& appl, const position& pos, std::size_t index)
{
  // These two conditions check pos' <= pos.
  if (pos.empty() || index >= pos.size())
  {
    return false;
  }
  else if (is_variable(appl[pos[index] - 1]))
  {
    return true;
  }
  else
  {
    assert(index < appl.size());
    return has_variable_higher_impl(static_cast<const application&>(appl[index]), pos, index + 1);
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
  for (auto& old_equation : equations)
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
  std::size_t current_state = m_automaton.root();

  // The empty position is the first index.
  m_subterms[0] = term;

  while (true)
  {
    // Retrieve the subterm.
    const apma_state& state = m_automaton.label(current_state);
    data_expression subterm = m_subterms[state.position];

    if (PrintMatchSteps)
    {
      mCRL2log(info) << "Matching subterm " << subterm << "\n";
    }

    // sigma := sigma \cup { (x, a) | x in L(s0) }
    if (m_automaton.label(current_state).variable.defined())
    {
      m_matching_sigma[m_automaton.label(current_state).variable] = subterm;
    }

    bool found_transition = false;
    if (is_application(subterm))
    {
      const auto& appl = static_cast<const application&>(subterm);

      // If delta(s0, a) is defined for some term a followed by suffix t'.
      auto [s_prime, found] = m_automaton.transition(current_state, get_head_index(appl.head()));
      if (found)
      {
        if (PrintMatchSteps)
        {
          mCRL2log(info) << "Took transition from " << current_state << " to " << s_prime << " with label " << appl.head() << "\n";
        }

        found_transition = true;
        current_state = s_prime;

        // Insert the subterms onto the position mapping.
        auto it = m_automaton.label(current_state).argument_positions.begin();
        for (const atermpp::aterm& argument : appl)
        {
          assert(it != m_automaton.label(current_state).argument_positions.end());
          m_subterms[*it] = static_cast<data_expression>(argument);

          ++it;
        }
      }
    }
    else if (is_function_symbol(subterm))
    {
      // If delta(s0, a) is defined for some term a followed by suffix t'.
      auto [s_prime, found] = m_automaton.transition(current_state, get_head_index(subterm));
      if (found)
      {
        if (PrintMatchSteps)
        {
          mCRL2log(info) << "Took transition " << current_state << " to " << s_prime << " with label " << static_cast<function_symbol>(subterm) << "\n";
        }

        found_transition = true;
        current_state = s_prime;
      }
    }

    if (!found_transition)
    {
      // If delta(s0, omega) is defined then
      auto [s_prime, found] = m_automaton.transition(current_state, m_not_equal_index);
      if (found)
      {
        // PMAMatch(delta(s0, omega)), t', sigma)
        if (PrintMatchSteps)
        {
          mCRL2log(info) << "Took transition from " << current_state << " to " << s_prime << " with label omega.\n";
        }

        current_state = s_prime;
      }
      // else return (emptyset, emptyset).
      else if (PrintMatchSteps)
      {
        mCRL2log(info) << "Matching failed, deadlock.\n";
        return;
      }
    }
  }

  if (PrintMatchSteps)
  {
    mCRL2log(info) << "Matching succeeded.\n";
  }

  m_match_set = &m_automaton.label(current_state).match_set;
}

template<typename Substitution>
const extended_data_equation* AdaptiveMatcher<Substitution>::next(Substitution& matching_sigma)
{
  if (m_match_set != nullptr)
  {
    while (m_match_index < m_match_set->size())
    {
      matching_sigma = m_matching_sigma;
      std::reference_wrapper<const linear_data_equation>& result = (*m_match_set)[m_match_index];
      ++m_match_index;

      if (!is_consistent(result.get().partition(), matching_sigma))
      {
        // This rule matched, but its variables are not consistent w.r.t. the substitution.
        continue;
      }

      return &result.get();
    }
  }

  return nullptr;
}

// Private functions

template<typename Substitution>
void AdaptiveMatcher<Substitution>::construct_apma(std::size_t s, data_expression pref)
{
  if (PrintConstructionSteps) { mCRL2log(info) << "pref = " << static_cast<atermpp::aterm>(pref) << "\n"; }

  // L := {l in L | l unifies with pref}
  std::vector<std::reference_wrapper<const linear_data_equation>> L;
  for (const auto& equation : m_linear_equations)
  {
    if (unify(equation.equation().lhs(), pref))
    {
      L.emplace_back(equation);
    }
  }

  // F := fringe(pref) \ intersection_{l in L} fringe(l) (the last part is ignored for now to ensure that the automaton is complete).
  std::set<position> F = fringe(pref);

  // if F = emptyset
  if (F.empty())
  {
    // M := M[L := L[s -> L]
    m_automaton.label(s).match_set.insert(m_automaton.label(s).match_set.begin(), L.begin(), L.end());
    return; // return M
  }
  else
  {
    // pos := select(F). For now the left-most position (left-to-right automaton).
    position pos = *std::min_element(F.begin(), F.end(), less_than);

    // M := M[L := L[s -> pos]]
    m_automaton.label(s).position = m_positions.insert(pos).first;

    // for f in F s.t. exist l in L : head(l[pos]) = f
    std::set<std::pair<mcrl2::data::function_symbol, std::size_t>> symbols;

    if (PrintConstructionSteps) { mCRL2log(info) << " Function symbols at selected position (" << pos << "):\n"; }
    for (const auto& equation : L)
    {
      // t := l[pos]
      std::optional<data_expression> t = at_position(static_cast<const application&>(equation.get().equation().lhs()), pos);

      if (t && (is_function_symbol(t.value()) || is_application_robust(t.value())))
      {
        // head(l[pos]) in F.
        symbols.insert(std::make_pair(static_cast<const function_symbol&>(get_head(t.value())), get_arity(t.value())));
      }
    }

    // This ensures that no duplicates are printed.
    if (PrintConstructionSteps)
    {
      for (const auto& [symbol, arity] : symbols)
      {
        mCRL2log(info) << "  " << symbol << "\n";
      }
    }

    for (const auto& [symbol, arity] : symbols)
    {
      // M := M[S := (S cup {s'})], where s' is a fresh state.
      std::size_t s_prime = m_automaton.add_state();

      // M := M[delta := (delta := delta(s, f) -> s')
      m_automaton.add_transition(s, mcrl2::core::index_traits<mcrl2::data::function_symbol, function_symbol_key_type, 2>::index(symbol), s_prime);

      // pref := f(omega^ar(f)), but we store in the correct position variables directly in the prefix.
      std::vector<data_expression> arguments;

      if (arity > 0)
      {
        // In this case there are multiple variables are introduced, their position is derived from current position.i where i is their
        // position in the application f(x_0, ..., x_n).
        position var_position = pos;
        var_position.push_back(0);
        for (std::size_t index = 0; index < arity; ++index)
        {
          var_position.back() = index + 1;
          arguments.push_back(position_variable(var_position));

          // Here, we also ensure that the arguments of this symbol are stored in the subterm table during evaluation.
          m_automaton.label(s_prime).argument_positions.emplace_back(m_positions.index(var_position));
        }

        application appl(symbol, arguments.begin(), arguments.end());

        construct_apma(s_prime, assign_at_position(pref, pos, appl));
      }
      else
      {
        // In this case the symbol is just a function symbol.
        construct_apma(s_prime, assign_at_position(pref, pos, symbol));
      }
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

      construct_apma(s_prime, assign_at_position(pref, pos, not_equal()));
    }
  }
}

// Explicit instantiations.

template class mcrl2::data::detail::AdaptiveMatcher<mutable_indexed_substitution<>>;
