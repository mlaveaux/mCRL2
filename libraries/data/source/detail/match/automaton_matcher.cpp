// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include "mcrl2/data/detail/match/automaton_matcher.h"

#include "mcrl2/data/detail/rewrite/jitty_jittyc.h"

#include <vector>

/// \brief Print the intermediate matches that succeeded.
constexpr bool PrintMatchSteps = false;

/// \brief Print the intermediate steps performed during construction.
constexpr bool PrintConstructionSteps = false;

using namespace mcrl2::data;
using namespace mcrl2::data::detail;

using namespace mcrl2::log;

void convert_to_string(const data_expression& expression, Pattern& pattern)
{
  if (is_application(expression))
  {
    const application& appl = static_cast<const application&>(expression);

    // Add the head symbol of appl.
    pattern.push_back(appl);

    for (const data_expression& argument : appl)
    {
      convert_to_string(argument, pattern);
    }
  }
  else
  {
    assert(is_function_symbol(expression) || is_variable(expression));
    pattern.push_back(expression);
  }
}

/// \returns A string representation (from sigma*) for a given term.
Pattern convert_to_string(const data_expression& expression)
{
  Pattern m_string;
  convert_to_string(expression, m_string);
  return m_string;
}

/// \returns The set omega / P, which are all patterns from P that start with a variable.
std::set<variable> get_head_variables(const PatternSet& P)
{
  std::set<variable> result;

  for (const Pattern& pattern : P)
  {
    if (is_variable(pattern.front()))
    {
      result.insert(static_cast<variable>(pattern.front()));
    }
  }

  return result;
}

/// \returns The set L / a, which are the tails of all patterns in P that start with a.
PatternSet a_derivatives(const PatternSet& P, const data_expression& a)
{
  PatternSet result;

  for (const Pattern& pattern : P)
  {
    const data_expression& head = pattern.front();
    const auto& appl = static_cast<const application&>(head);

    if (head == a || (is_application(head) && head.size() > 0 && appl.head() == a))
    {
      Pattern tail(pattern.begin() + 1, pattern.end());
      result.insert(tail);
    }
  }

  return result;
}

/// \returns The set L / omega, which is the union of all L / x for x in V (the set of all variables).
PatternSet omega_derivatives(const PatternSet& P)
{
  PatternSet result;

  for (const variable& variable : get_head_variables(P))
  {
    PatternSet P_x = a_derivatives(P, variable);
    result.insert(P_x.begin(), P_x.end());
  }

  return result;
}

// Public functions

template<typename Substitution>
AutomatonMatcher<Substitution>::AutomatonMatcher(const data_equation_vector& equations)
{
  // Compute the root by converting all terms to string representations.
  PatternSet root;

  for (auto& equation : equations)
  {
    Pattern string = convert_to_string(equation.lhs());

    // Add the index of the equation
    string.push_back(end_of_string(equation));
    m_mapping.insert(std::make_pair(equation,
      data_equation_extended(equation,
        ConstructionStack(equation.condition()),
        ConstructionStack(equation.rhs())))
    );

    root.insert(string);
  }

  // Construct the automaton.
  m_root_state = add_fresh_state();
  construct_rec(root, m_root_state);
}

template<typename Substitution>
void AutomatonMatcher<Substitution>::match(const data_expression& term)
{
  m_matching_sigma.clear();
  m_match_set = nullptr;
  m_match_index = 0;

  // Start with the root state.
  pma_state* current_state = m_root_state;

  // Implemented iteratively to return the matching set.
  std::stack<data_expression> m_stack;
  m_stack.push(term);

  while (!m_stack.empty())
  {
    data_expression subterm = m_stack.top();
    m_stack.pop();

    if (PrintMatchSteps)
    {
      mCRL2log(info) << "Matching subterm " << subterm << "\n";
    }

    // sigma := sigma \cup { (x, a) | x in L(s0) }
    for (const variable& variable : current_state->variables)
    {
      m_matching_sigma[variable] = term;
    }

    bool found_transition = false;
    if (is_application(subterm))
    {
      const auto& appl = static_cast<const application&>(subterm);

      // If delta(s0, a) is defined for some term a followed by suffix t'.
      auto result = m_transitions.find(std::make_pair(current_state, appl.head()));
      if (result != m_transitions.end())
      {
        if (PrintMatchSteps)
        {
          mCRL2log(info) << "Took transition from " << current_state << " to " << (*result).second << " with label " << appl.head() << "\n";
        }

        found_transition = true;
        current_state = (*result).second;

        // Insert in reverse order.
        for (int index = static_cast<int>(appl.size() - 1); index >= 0; --index)
        {
          assert(index >= 0);
          m_stack.push(appl[static_cast<std::size_t>(index)]);
        }
      }
    }
    else if (is_function_symbol(subterm))
    {
      const auto& symbol = static_cast<const function_symbol&>(subterm);

      // If delta(s0, a) is defined for some term a followed by suffix t'.
      auto result = m_transitions.find(std::make_pair(current_state, symbol));
      if (result != m_transitions.end())
      {
        if (PrintMatchSteps)
        {
          mCRL2log(info) << "Took transition " << current_state << " to " << (*result).second << " with label " << symbol << "\n";
        }

        found_transition = true;
        current_state = (*result).second;
      }
    }

    if (!found_transition)
    {
      // If delta(s0, omega) is defined then
      auto result = m_transitions.find(std::make_pair(current_state, m_omega));
      if (result != m_transitions.end())
      {
        // PMAMatch(delta(s0, omega)), t', sigma)
        if (PrintMatchSteps)
        {
          mCRL2log(info) << "Took transition from " << current_state << " to " << (*result).second << " with omega.\n";
        }

        found_transition = true;
        current_state = (*result).second;
      }
    }

    if (!found_transition)
    {
      // else return (emptyset, emptyset).
      if (PrintMatchSteps)
      {
        mCRL2log(info) << "Matching failed, deadlock.\n";
      }
    }
  }

  if (PrintMatchSteps)
  {
    mCRL2log(info) << "Matching succeeded.\n";
  }

  m_match_set = &current_state->match_set;
}

template<typename Substitution>
const data_equation_extended* AutomatonMatcher<Substitution>::next(Substitution& matching_sigma)
{
  matching_sigma = m_matching_sigma;

  if (m_match_set != nullptr)
  {
    if (m_match_index < m_match_set->size())
    {
      return m_match_set[m_match_index++];
    }
  }

  return nullptr;
}

// Private functions

template<typename Substitution>
typename AutomatonMatcher<Substitution>::pma_state* AutomatonMatcher<Substitution>::add_fresh_state()
{
  // Create a fresh state s0
  m_states.push_back(std::unique_ptr<pma_state>(new pma_state()));
  return m_states.back().get();
}

template<typename Substitution>
void AutomatonMatcher<Substitution>::construct_rec(const PatternSet& L, pma_state* s0)
{
  if (PrintConstructionSteps)
  {
    mCRL2log(info) << "L = ";
    print_pattern_set(L);
  }

  // L := (s0, omega / L) can be represented by giving s0 the label omega / L.
  s0->variables = get_head_variables(L);

  // If L = { e_i | i in I for some I subset N }, i.e. it only consists of final matches.
  std::vector<std::reference_wrapper<const data_equation_extended>> match_set;

  for (const Pattern& p : L)
  {
    if (is_end_of_string(p.front()))
    {
      assert(p.size() == 1); // Should only consist of the empty string.

      auto result = m_mapping.find(static_cast<const end_of_string&>(p.front()).equation());
      assert(result != m_mapping.end());
      match_set.push_back((*result).second);
    }
  }

  if (!match_set.empty())
  {
    // Return F[s0] := {r in R | the epsilon of r is labelled with some i in I }
    s0->match_set = match_set;
    return;
  }

  // L_omega = L / omega.
  PatternSet L_omega = omega_derivatives(L);

  if (PrintConstructionSteps)
  {
    mCRL2log(info) << "L / omega = ";
    print_pattern_set(L_omega);
  }

  // For f in Sigma \ V such that L / f != empty.
  for (const data_expression& f : get_head_symbols(L))
  {
    assert(!is_end_of_string(f));
    assert(!is_omega(f));

    if (is_application(f))
    {
      const auto& appl = static_cast<const application&>(f);

      // L_f = L / f union { omega^ar(f) p | p in L_omega }
      PatternSet L_f = a_derivatives(L, appl.head());

      PatternSet L2 = prepend_omegas(L_omega, appl);

      L_f.insert(L2.begin(), L2.end());

      // RecConstruct(L_f, R)
      pma_state* s0_prime = add_fresh_state();

      if (PrintConstructionSteps)
      {
        mCRL2log(info) << "L / f cup { omega^ar(f) p | p in L / omega } with f equal to " << appl.head() << " = ";
        print_pattern_set(L_f);
        mCRL2log(info) << "Added transition from " << s0 << " to " << s0_prime << " with label " << appl.head() << ".\n";
      }

      construct_rec(L_f, s0_prime);

      // delta := delta cup delta' cup { (s0, f, s'0) }
      m_transitions.insert(std::make_pair(std::make_pair(s0, appl.head()), s0_prime));
    }
    else
    {
      assert(is_function_symbol(f));
      PatternSet L_f = a_derivatives(L, f);

      // RecConstruct(L_f, R)
      pma_state* s0_prime = add_fresh_state();

      if (PrintConstructionSteps)
      {
        mCRL2log(info) << "L / " << f << " = ";
        print_pattern_set(L_f);
        mCRL2log(info) << "Added transition from " << s0 << " to " << s0_prime << " with label " << f << ".\n";
      }

      construct_rec(L_f, s0_prime);

      // delta := delta cup delta' cup { (s0, f, s'0) }
      m_transitions.insert(std::make_pair(std::make_pair(s0, f), s0_prime));
    }

    // the automata are merged by having a global transition system and a labelling local to each state.
  }

  if (!L_omega.empty())
  {
    // RecConstruct(L_omega, R)
    pma_state* s0_prime = add_fresh_state();
    construct_rec(L_omega, s0_prime);

    // delta := delta cup delta' cup { (s0, omega, s'0) }
    if (PrintConstructionSteps)
    {
      mCRL2log(info) << "Added transition from " << s0 << " to " << s0_prime << " with label omega.\n";
    }

    m_transitions.insert(std::make_pair(std::make_pair(s0, m_omega), s0_prime));

    // the automata are merged by having a global transition system and a labelling local to each state.
  }

  // The automaton is not returned, but constructed in place.
  if (PrintConstructionSteps)
  {
    mCRL2log(info) << "Returning.\n";
  }
}

template<typename Substitution>
std::set<data_expression> AutomatonMatcher<Substitution>::get_head_symbols(const PatternSet& P) const
{
  std::set<data_expression> result;

  for (const Pattern& pattern : P)
  {
    if (!is_variable(pattern.front()) && !is_omega(pattern.front()))
    {
      result.insert(pattern.front());
    }
  }

  return result;
}

template<typename Substitution>
PatternSet AutomatonMatcher<Substitution>::prepend_omegas(const PatternSet& L, const application& f) const
{
  PatternSet result;

  for (Pattern pattern : L)
  {
    // Prepend arity(f) omega symbols.
    for (std::size_t i = 0; i < f.size(); ++i)
    {
      pattern.emplace(pattern.begin(), omega());
    }

    result.insert(pattern);
  }

  return result;
}

template<typename Substitution>
void AutomatonMatcher<Substitution>::print_pattern_set(const PatternSet& set) const
{
  mCRL2log(info) << "{ ";
  bool first = true;
  for (const auto& pattern : set)
  {
    if (!first)
    {
      mCRL2log(info) << ", ";
    }

    for (const auto& expression : pattern)
    {
      if (is_omega(expression))
      {
        mCRL2log(info) << " omega";
      }
      else if (!is_end_of_string(expression))
      {
        mCRL2log(info) << " " << expression;
      }
    }
    first = false;
  }

  mCRL2log(info) << "} \n";
}
