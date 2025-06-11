// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pbes/pbes.h
/// \brief The class pbes.

#ifndef MCRL2_PBES_PBES_H
#define MCRL2_PBES_PBES_H

#include "mcrl2/data/detail/equal_sorts.h"
#include "mcrl2/pbes/pbes_equation.h"

namespace mcrl2
{

/// \brief The main namespace for the PBES library.
namespace pbes_system
{

class pbes;
void complete_data_specification(pbes&);

// template function overloads
std::string pp(const pbes& x, bool precedence_aware = true);
void normalize_sorts(pbes& x, const data::sort_specification& sortspec);
void translate_user_notation(pbes_system::pbes& x);
std::set<data::sort_expression> find_sort_expressions(const pbes_system::pbes& x);
std::set<data::variable> find_all_variables(const pbes_system::pbes& x);
std::set<data::variable> find_free_variables(const pbes_system::pbes& x);
std::set<data::function_symbol> find_function_symbols(const pbes_system::pbes& x);

bool is_well_typed_equation(const pbes_equation& eqn,
                            const std::set<data::sort_expression>& declared_sorts,
                            const std::set<data::variable>& declared_global_variables,
                            const data::data_specification& data_spec
                           );

bool is_well_typed_pbes(const std::set<data::sort_expression>& declared_sorts,
                        const std::set<data::variable>& declared_global_variables,
                        const std::set<data::variable>& occurring_global_variables,
                        const std::set<propositional_variable>& declared_variables,
                        const std::set<propositional_variable_instantiation>& occ,
                        const propositional_variable_instantiation& init,
                        const data::data_specification& data_spec
                       );

atermpp::aterm pbes_to_aterm(const pbes& p);

/// \brief parameterized boolean equation system
// <PBES>         ::= PBES(<DataSpec>, <GlobVarSpec>, <PBEqnSpec>, <PBInit>)
// <PBEqnSpec>    ::= PBEqnSpec(<PBEqn>*)
class pbes
{
  public:
    typedef pbes_equation equation_type;

  protected:
    /// \brief The data specification
    data::data_specification m_data;

    /// \brief The sequence of pbes equations
    std::vector<pbes_equation> m_equations;

    /// \brief The set of global variables
    std::set<data::variable> m_global_variables;

    /// \brief The initial state
    propositional_variable_instantiation m_initial_state;

    /// \brief Returns the predicate variables appearing in the left hand side of an equation.
    /// \return The predicate variables appearing in the left hand side of an equation.
    std::set<propositional_variable> compute_declared_variables() const
    {
      std::set<propositional_variable> result;
      for (const pbes_equation& eqn: equations())
      {
        result.insert(eqn.variable());
      }
      return result;
    }

    /// \brief Checks if the propositional variable instantiation v appears with the right type in the
    /// sequence of propositional variable declarations [first, last).
    /// \param first Start of a sequence of propositional variable declarations
    /// \param last End of a sequence of propositional variable declarations
    /// \param v A propositional variable instantation
    /// \param data_spec A data specification.
    /// \return True if the type of \p v is matched correctly
    /// \param v A propositional variable instantiation
    template <typename Iter>
    bool is_declared_in(Iter first, Iter last, const propositional_variable_instantiation& v, const data::data_specification& data_spec) const
    {
      for (Iter i = first; i != last; ++i)
      {
        if (i->name() == v.name() && data::detail::equal_sorts(i->parameters(), v.parameters(), data_spec))
        {
          return true;
        }
      }
      return false;
    }

  public:
    /// \brief Constructor.
    pbes() = default;

    /// \brief Constructor.
    /// \param data A data specification
    /// \param equations A sequence of pbes equations
    /// \param initial_state A propositional variable instantiation
    pbes(const data::data_specification& data,
         const std::vector<pbes_equation>& equations,
         propositional_variable_instantiation initial_state)
      :
      m_data(data),
      m_equations(equations),
      m_initial_state(initial_state)
    {
      m_global_variables = pbes_system::find_free_variables(*this);
      assert(core::detail::check_rule_PBES(pbes_to_aterm(*this)));
      assert(is_well_typed());
    }

    /// \brief Constructor.
    /// \param data A data specification
    /// \param equations A sequence of pbes equations
    /// \param global_variables A sequence of free variables
    /// \param initial_state A propositional variable instantiation
    pbes(const data::data_specification& data,
         const std::set<data::variable>& global_variables,
         const std::vector<pbes_equation>& equations,
         propositional_variable_instantiation initial_state)
      :
      m_data(data),
      m_equations(equations),
      m_global_variables(global_variables),
      m_initial_state(initial_state)
    {
      assert(core::detail::check_rule_PBES(pbes_to_aterm(*this)));
      assert(is_well_typed());
    }

    /// \brief Returns the data specification.
    /// \return The data specification of the pbes
    const data::data_specification& data() const
    {
      return m_data;
    }

    /// \brief Returns the data specification.
    /// \return The data specification of the pbes
    data::data_specification& data()
    {
      return m_data;
    }

    /// \brief Returns the equations.
    /// \return The equations.
    const std::vector<pbes_equation>& equations() const
    {
      return m_equations;
    }

    /// \brief Returns the equations.
    /// \return The equations.
    std::vector<pbes_equation>& equations()
    {
      return m_equations;
    }

    /// \brief Returns the declared free variables of the pbes.
    /// \return The declared free variables of the pbes.
    const std::set<data::variable>& global_variables() const
    {
      return m_global_variables;
    }

    /// \brief Returns the declared free variables of the pbes.
    /// \return The declared free variables of the pbes.
    std::set<data::variable>& global_variables()
    {
      return m_global_variables;
    }

    /// \brief Returns the initial state.
    /// \return The initial state.
    const propositional_variable_instantiation& initial_state() const
    {
      return m_initial_state;
    }

    /// \brief Returns the initial state.
    /// \return The initial state.
    propositional_variable_instantiation& initial_state()
    {
      return m_initial_state;
    }

    /// \brief Returns the set of binding variables of the pbes.
    /// This is the set variables that occur on the left hand side of an equation.
    /// \return The set of binding variables of the pbes.
    std::set<propositional_variable> binding_variables() const
    {
      using namespace std::rel_ops; // for definition of operator!= in terms of operator==

      std::set<propositional_variable> result;
      for (const pbes_equation& eqn: equations())
      {
        result.insert(eqn.variable());
      }
      return result;
    }

    /// \brief Returns the set of occurring propositional variable instantiations of the pbes.
    /// This is the set of variables that occur in the right hand side of an equation.
    /// \return The occurring propositional variable instantiations of the pbes
    std::set<propositional_variable_instantiation> occurring_variable_instantiations() const;

    /// \brief Returns the set of occurring propositional variable declarations of the pbes, i.e.
    /// the propositional variable declarations that occur in the right hand side of an equation.
    /// \return The occurring propositional variable declarations of the pbes
    std::set<propositional_variable> occurring_variables() const
    {
      std::set<propositional_variable> result;
      std::set<propositional_variable_instantiation> occ = occurring_variable_instantiations();
      std::map<core::identifier_string, propositional_variable> declared_variables;
      for (const pbes_equation& eqn: equations())
      {
        declared_variables[eqn.variable().name()] = eqn.variable();
      }
      for (const propositional_variable_instantiation& v: occ)
      {
        result.insert(declared_variables[v.name()]);
      }
      return result;
    }

    /// \brief True if the pbes is closed
    /// \return Returns true if all occurring variables are binding variables, and the initial state variable is a binding variable.
    bool is_closed() const
    {
      std::set<propositional_variable> bnd = binding_variables();
      std::set<propositional_variable> occ = occurring_variables();
      return std::includes(bnd.begin(), bnd.end(), occ.begin(), occ.end()) && is_declared_in(bnd.begin(), bnd.end(), initial_state(), data());
    }

    /// \brief Checks if the PBES is well typed
    /// \return True if
    /// <ul>
    /// <li>the sorts occurring in the free variables of the equations are declared in the data specification</li>
    /// <li>the sorts occurring in the binding variable parameters are declared in the data specification </li>
    /// <li>the sorts occurring in the quantifier variables of the equations are declared in the data specification </li>
    /// <li>the binding variables of the equations have unique names (well formedness)</li>
    /// <li>the global variables occurring in the equations are declared in global_variables()</li>
    /// <li>the global variables occurring in the equations with the same name are identical</li>
    /// <li>the declared global variables and the quantifier variables occurring in the equations have different names</li>
    /// <li>the predicate variable instantiations occurring in the equations match with their declarations</li>
    /// <li>the predicate variable instantiation occurring in the initial state matches with the declaration</li>
    /// <li>the data specification is well typed</li>
    /// </ul>
    /// N.B. Conflicts between the types of instantiations and declarations of binding variables are not checked!
    bool is_well_typed() const
    {
      std::set<data::sort_expression> declared_sorts = data::detail::make_set(data().sorts());
      const std::set<data::variable>& declared_global_variables = global_variables();
      std::set<data::variable> occurring_global_variables = pbes_system::find_free_variables(*this);
      std::set<propositional_variable> declared_variables = compute_declared_variables();
      std::set<propositional_variable_instantiation> occ = occurring_variable_instantiations();

      // check 1), 4), 5), 6), 8) and 9)
      if (!is_well_typed_pbes(declared_sorts, declared_global_variables, occurring_global_variables, declared_variables, occ, initial_state(), data()))
      {
        return false;
      }

      // check 2), 3) and 7)
      for (const pbes_equation& eqn: equations())
      {
        if (!is_well_typed_equation(eqn, declared_sorts, declared_global_variables, data()))
        {
          return false;
        }
      }

      // check 10)
      return data().is_well_typed();
    }
};

//--- start generated class pbes ---//
// prototype declaration
std::string pp(const pbes& x, bool precedence_aware);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const pbes& x)
{
  return out << pbes_system::pp(x);
}
//--- end generated class pbes ---//


/// \brief Adds all sorts that appear in the PBES \a p to the data specification of \a p.
/// \param p a PBES.
inline
void complete_data_specification(pbes& p)
{
  std::set<data::sort_expression> s = pbes_system::find_sort_expressions(p);
  p.data().add_context_sorts(s);
}

/// \brief Equality operator on PBESs
/// \return True if the PBESs have exactly the same internal representation. Note
/// that this is in general not a very useful test.
// TODO: improve the comparison
inline
bool operator==(const pbes& p1, const pbes& p2)
{
  return pbes_to_aterm(p1) == pbes_to_aterm(p2);
}

} // namespace pbes_system

} // namespace mcrl2

#endif // MCRL2_PBES_PBES_H
