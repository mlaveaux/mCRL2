// Author(s): Jan Friso Groote. Based on pbes_equation.h by Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pres/pres_equation.h
/// \brief The class pres_equation.

#ifndef MCRL2_PRES_PRES_EQUATION_H
#define MCRL2_PRES_PRES_EQUATION_H

#include "mcrl2/data/detail/sequence_algorithm.h"
#include "mcrl2/pbes/fixpoint_symbol.h"
#include "mcrl2/pres/pres_expression.h"

namespace mcrl2
{

namespace pres_system
{

typedef pbes_system::fixpoint_symbol fixpoint_symbol;

class pres_equation;
atermpp::aterm pres_equation_to_aterm(const pres_equation& eqn);
bool is_well_typed(const pres_equation& eqn);
bool has_propositional_variables(const pres_expression& x);

/// \brief pres equation.
class pres_equation
{
  protected:
    /// \brief The fixpoint symbol of the equation
    fixpoint_symbol m_symbol;

    /// \brief The variable on the left hand side of the equation
    propositional_variable m_variable;

    /// \brief The expression on the right hand side of the equation
    pres_expression m_formula;

  public:
    /// \brief The expression type of the equation.
    typedef pres_expression term_type;

    /// \brief The variable type of the equation.
    typedef propositional_variable variable_type;

    /// \brief The symbol type of the equation.
    typedef fixpoint_symbol symbol_type;

    /// \brief Constructor.
    pres_equation() = default;

    /// \brief Constructor.
    /// \param symbol A fixpoint symbol
    /// \param variable A propositional variable declaration
    /// \param expr A PRES expression
    pres_equation(const fixpoint_symbol& symbol, const propositional_variable& variable, const pres_expression& expr)
      :
      m_symbol(symbol),
      m_variable(variable),
      m_formula(expr)
    {
    }

    /// \brief Returns the fixpoint symbol of the equation.
    /// \return The fixpoint symbol of the equation.
    const fixpoint_symbol& symbol() const
    {
      return m_symbol;
    }

    /// \brief Returns the fixpoint symbol of the equation.
    /// \return The fixpoint symbol of the equation.
    fixpoint_symbol& symbol()
    {
      return m_symbol;
    }

    /// \brief Returns the pres variable of the equation.
    /// \return The pres variable of the equation.
    const propositional_variable& variable() const
    {
      return m_variable;
    }

    /// \brief Returns the pres variable of the equation.
    /// \return The pres variable of the equation.
    propositional_variable& variable()
    {
      return m_variable;
    }

    /// \brief Returns the predicate formula on the right hand side of the equation.
    /// \return The predicate formula on the right hand side of the equation.
    const pres_expression& formula() const
    {
      return m_formula;
    }

    /// \brief Returns the predicate formula on the right hand side of the equation.
    /// \return The predicate formula on the right hand side of the equation.
    pres_expression& formula()
    {
      return m_formula;
    }

    /// \brief Returns true if the predicate formula on the right hand side contains no predicate variables.
    // (Comment Wieger: is_const would be a better name)
    /// \return True if the predicate formula on the right hand side contains no predicate variables.
    bool is_solved() const;

    /// \brief Swaps the contents
    void swap(pres_equation& other)
    {
      using std::swap;
      swap(m_symbol, other.m_symbol);
      swap(m_variable, other.m_variable);
      swap(m_formula, other.m_formula);
    }
   
    void mark(std::stack<std::reference_wrapper<atermpp::detail::_aterm>>& todo) const
    {
      mark_term(m_symbol, todo);
    }
};

//--- start generated class pres_equation ---//
/// \\brief list of pres_equations
typedef atermpp::term_list<pres_equation> pres_equation_list;

/// \\brief vector of pres_equations
typedef std::vector<pres_equation>    pres_equation_vector;

// prototype declaration
std::string pp(const pres_equation& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const pres_equation& x)
{
  return out << pres_system::pp(x);
}

/// \\brief swap overload
inline void swap(pres_equation& t1, pres_equation& t2)
{
  t1.swap(t2);
}
//--- end generated class pres_equation ---//

inline bool
operator==(const pres_equation& x, const pres_equation& y)
{
  return x.symbol() == y.symbol() &&
         x.variable() == y.variable() &&
         x.formula() == y.formula();
}

inline bool
operator!=(const pres_equation& x, const pres_equation& y)
{
  return !(x == y);
}

/// \brief Conversion to atermaPpl.
/// \return The specification converted to aterm format.
inline
atermpp::aterm pres_equation_to_aterm(const pres_equation& eqn)
{
  return atermpp::aterm(core::detail::function_symbol_PREqn(), eqn.symbol(), eqn.variable(), eqn.formula());
}

// Overload for pp for propositional variables.
std::string pp(const propositional_variable& v);

// template function overloads
std::string pp(const pres_equation_vector& x);
void normalize_sorts(pres_equation_vector& x, const data::sort_specification& sortspec);
std::set<data::variable> find_free_variables(const pres_system::pres_equation& x);

} // namespace pres_system

} // namespace mcrl2

#endif // MCRL2_PRES_PRES_EQUATION_H
