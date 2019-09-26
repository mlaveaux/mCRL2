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
#include "mcrl2/utilities/stopwatch.h"
#include "mcrl2/core/index_traits.h"

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

/// \returns A unique index for the head symbol that the given term starts with.
inline std::size_t get_head_index(const data_expression& term)
{
  return mcrl2::core::index_traits<mcrl2::data::function_symbol, function_symbol_key_type, 2>::index(static_cast<const function_symbol&>(get_nested_head(term)));
}

using position = std::vector<std::size_t>;


/// \returns A position [1,2,3,4] as the string 1.2.3.4
std::string to_text(const position& position)
{
  std::string result;

  bool first = true;
  for (std::size_t index : position)
  {
    if (!first)
    {
      result += ".";
    }

    result += std::to_string(index);

    first = false;
  }

  return result;
}

/// \brief Renames every variable to a unique name by using its position as identifier.
template<typename Substitution>
void rename_variables_position(const atermpp::aterm_appl& appl, position current, Substitution& sigma)
{
  if (is_variable(appl))
  {
    const auto& var = variable{appl};
    sigma[var] = variable(mcrl2::core::identifier_string(to_text(current)), var.sort());
  }
  else
  {
    current.emplace_back(0);
    for (const atermpp::aterm& argument : appl)
    {
      rename_variables_position(static_cast<const atermpp::aterm_appl&>(argument), current, sigma);
      ++current.back();
    }
  }
}

template<typename Substitution>
void rename_variables_position(const atermpp::aterm_appl& appl, Substitution& sigma)
{
  rename_variables_position(appl, position(), sigma);
}

/// \brief Rename the variables in the data_equation such that at each position they are unique and adapt the equivalence classes.
std::pair<data_equation, partition> rename_variables_unique(std::pair<data_equation, partition> rules)
{
  mutable_indexed_substitution<variable, variable> sigma;
  rename_variables_position(rules.first.lhs(), sigma);

  // Rename all variables in the resulting partition to the name indicated by sigma.
  partition result;
  for (const std::vector<variable>& set : rules.second)
  {
    std::vector<variable> new_vars;
    for (const variable& var : set)
    {
      new_vars.push_back(sigma(var));
    }
    result.push_back(new_vars);
  }

  return std::make_pair(replace_variables(rules.first, sigma), result);
}

// Public functions

template<typename Substitution>
AutomatonMatcher<Substitution>::AutomatonMatcher(const data_equation_vector& equations)
  : m_automaton()
{
  mcrl2::utilities::stopwatch construction;
  enumerator_identifier_generator generator("@");

  // Compute the root by converting all terms to string representations.
  PatternSet root;

  for (auto& old_equation : equations)
  {
    // Rename the variables in the equation
    std::pair<data_equation, partition> result = rename_variables_unique(make_linear(old_equation, generator));

    auto& equation = result.first;
    if (!result.second.empty())
    {
      // For now, print the now linear equations.
      mCRL2log(info) << "Renamed non-linear equation " << equation << ".\n";
    }

    Pattern string = convert_to_string(equation.lhs());

    // Add the index of the equation
    string.push_back(end_of_string(equation));
    m_mapping.insert(std::make_pair(equation, std::make_pair(data_equation_extended(equation), result.second)));

    root.insert(string);
  }

  // Determine the index of omega.
  m_omega_index = mcrl2::core::index_traits<mcrl2::data::function_symbol, function_symbol_key_type, 2>::index(static_cast<const function_symbol&>(omega()));

  // Construct the automaton.
  construct_rec(root, m_automaton.root());

  mCRL2log(info) << "Matching automaton (states: " << m_automaton.states() << ", transitions: " << m_automaton.transitions() << ") construction took " << construction.time() << " milliseconds.\n";
}

template<typename Substitution>
void AutomatonMatcher<Substitution>::match(const data_expression& term)
{
  m_matching_sigma.clear();
  m_match_set = nullptr;
  m_match_index = 0;

  // Start with the root state.
  std::size_t current_state = m_automaton.root();

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
    if (m_automaton.label(current_state).variable.defined())
    {
      m_matching_sigma[m_automaton.label(current_state).variable] = subterm;
    }

    bool found_transition = false;
    if (is_application(subterm))
    {
      const auto& appl = static_cast<const application&>(subterm);

      // If delta(s0, a) is defined for some term a followed by suffix t'.
      auto result = m_automaton.transition(current_state, get_head_index(appl.head()));
      if (result.second)
      {
        if (PrintMatchSteps)
        {
          mCRL2log(info) << "Took transition from " << current_state << " to " << result.first << " with label " << appl.head() << "\n";
        }

        found_transition = true;
        current_state = result.first;

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
      // If delta(s0, a) is defined for some term a followed by suffix t'.
      auto result = m_automaton.transition(current_state, get_head_index(subterm));
      if (result.second)
      {
        if (PrintMatchSteps)
        {
          mCRL2log(info) << "Took transition " << current_state << " to " << result.first << " with label " << static_cast<function_symbol>(subterm) << "\n";
        }

        found_transition = true;
        current_state = result.first;
      }
    }

    if (!found_transition)
    {
      // If delta(s0, omega) is defined then
      auto result = m_automaton.transition(current_state, m_omega_index);
      if (result.second)
      {
        // PMAMatch(delta(s0, omega)), t', sigma)
        if (PrintMatchSteps)
        {
          mCRL2log(info) << "Took transition from " << current_state << " to " << result.first << " with label omega.\n";
        }

        current_state = result.first;
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
const data_equation_extended* AutomatonMatcher<Substitution>::next(Substitution& matching_sigma)
{
  if (m_match_set != nullptr)
  {
    while (m_match_index < m_match_set->size())
    {
      matching_sigma = m_matching_sigma;
      std::reference_wrapper<const linear_data_equation>& result = (*m_match_set)[m_match_index];
      ++m_match_index;

      if (!is_consistent(result.get().second, matching_sigma))
      {
        // This rule matched, but its variables are not consistent w.r.t. the substitution.
        continue;
      }

      return &result.get().first;
    }
  }

  return nullptr;
}

// Private functions

template<typename Substitution>
void AutomatonMatcher<Substitution>::construct_rec(const PatternSet& L, std::size_t s0)
{
  if (PrintConstructionSteps)
  {
    mCRL2log(info) << "L = ";
    print_pattern_set(L);
  }

  // L := (s0, omega / L) can be represented by giving s0 the label omega / L.
  std::set<variable> variables = get_head_variables(L);

  if (variables.size() > 0)
  {
    assert(variables.size() == 1);
    m_automaton.label(s0).variable = *variables.begin();
  }

  // If L = { e_i | i in I for some I subset N }, i.e. it only consists of final matches.
  std::vector<std::reference_wrapper<const linear_data_equation>> match_set;

  for (const Pattern& p : L)
  {
    if (is_end_of_string(p.front()))
    {
      assert(p.size() == 1); // Should only consist of the empty string.

      auto result = m_mapping.find(static_cast<const end_of_string&>(p.front()).equation());
      assert(result != m_mapping.end());
      match_set.push_back(result->second);
    }
  }

  if (!match_set.empty())
  {
    // Return F[s0] := {r in R | the epsilon of r is labelled with some i in I }
    m_automaton.label(s0).match_set = match_set;
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
      std::size_t s0_prime = m_automaton.add_state();

      if (PrintConstructionSteps)
      {
        mCRL2log(info) << "L / f cup { omega^ar(f) p | p in L / omega } with f equal to " << appl.head() << " = ";
        print_pattern_set(L_f);
        mCRL2log(info) << "Added transition from " << s0 << " to " << s0_prime << " with label " << appl.head() << ".\n";
      }

      construct_rec(L_f, s0_prime);

      // delta := delta cup delta' cup { (s0, f, s'0) }
      m_automaton.add_transition(s0, get_head_index(appl.head()), s0_prime);
    }
    else
    {
      assert(is_function_symbol(f));
      PatternSet L_f = a_derivatives(L, f);

      // RecConstruct(L_f, R)
      std::size_t s0_prime = m_automaton.add_state();

      if (PrintConstructionSteps)
      {
        mCRL2log(info) << "L / " << f << " = ";
        print_pattern_set(L_f);
        mCRL2log(info) << "Added transition from " << s0 << " to " << s0_prime << " with label " << f << ".\n";
      }

      construct_rec(L_f, s0_prime);

      // delta := delta cup delta' cup { (s0, f, s'0) }
      m_automaton.add_transition(s0, get_head_index(f), s0_prime);
    }

    // the automata are merged by having a global transition system and a labelling local to each state.
  }

  if (!L_omega.empty())
  {
    // RecConstruct(L_omega, R)
    std::size_t s0_prime = m_automaton.add_state();
    construct_rec(L_omega, s0_prime);

    // delta := delta cup delta' cup { (s0, omega, s'0) }
    if (PrintConstructionSteps)
    {
      mCRL2log(info) << "Added transition from " << s0 << " to " << s0_prime << " with label omega.\n";
    }

    m_automaton.add_transition(s0, m_omega_index, s0_prime);

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

// Explicit instantiations.

template class mcrl2::data::detail::AutomatonMatcher<mutable_indexed_substitution<>>;
