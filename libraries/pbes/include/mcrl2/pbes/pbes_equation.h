// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pbes/pbes_equation.h
/// \brief The class pbes_equation.

#ifndef MCRL2_PBES_PBES_EQUATION_H
#define MCRL2_PBES_PBES_EQUATION_H

#include "mcrl2/data/detail/sequence_algorithm.h"
#include "mcrl2/pbes/fixpoint_symbol.h"
#include "mcrl2/pbes/pbes_expression.h"

namespace mcrl2
{

namespace pbes_system
{

class pbes_equation;
atermpp::aterm pbes_equation_to_aterm(const pbes_equation& eqn);
bool is_well_typed(const pbes_equation& eqn);
bool has_propositional_variables(const pbes_expression& x);

/// \brief pbes equation.
class pbes_equation
{
  protected:
    /// \brief The fixpoint symbol of the equation
    fixpoint_symbol m_symbol;

    /// \brief The variable on the left hand side of the equation
    propositional_variable m_variable;

    /// \brief The expression on the right hand side of the equation
    pbes_expression m_formula;

  public:
    /// \brief The expression type of the equation.
    typedef pbes_expression term_type;

    /// \brief The variable type of the equation.
    typedef propositional_variable variable_type;

    /// \brief The symbol type of the equation.
    typedef fixpoint_symbol symbol_type;

    /// \brief Constructor.
    pbes_equation() = default;

    /// \brief Constructor.
    /// \param symbol A fixpoint symbol
    /// \param variable A propositional variable declaration
    /// \param expr A PBES expression
    pbes_equation(const fixpoint_symbol& symbol, const propositional_variable& variable, const pbes_expression& expr)
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

    /// \brief Returns the pbes variable of the equation.
    /// \return The pbes variable of the equation.
    const propositional_variable& variable() const
    {
      return m_variable;
    }

    /// \brief Returns the pbes variable of the equation.
    /// \return The pbes variable of the equation.
    propositional_variable& variable()
    {
      return m_variable;
    }

    /// \brief Returns the predicate formula on the right hand side of the equation.
    /// \return The predicate formula on the right hand side of the equation.
    const pbes_expression& formula() const
    {
      return m_formula;
    }

    /// \brief Returns the predicate formula on the right hand side of the equation.
    /// \return The predicate formula on the right hand side of the equation.
    pbes_expression& formula()
    {
      return m_formula;
    }

    /// \brief Returns true if the predicate formula on the right hand side contains no predicate variables.
    // (Comment Wieger: is_const would be a better name)
    /// \return True if the predicate formula on the right hand side contains no predicate variables.
    bool is_solved() const;

    /// \brief Swaps the contents
    void swap(pbes_equation& other)
    {
      using std::swap;
      swap(m_symbol, other.m_symbol);
      swap(m_variable, other.m_variable);
      swap(m_formula, other.m_formula);
    }

    /// \brief A comparison operator on pbes equations.
    /// \detail The comparison is on the addresses of aterm objects and can therefore differ in every run.
    /// \parameter other The variable to compare with. 
    /// \return True if the this argument is smaller than other.
    bool operator<(const pbes_equation& other) const
    {
      return (m_formula < other.m_formula) ||
             ((m_formula == other.m_formula &&
                 (m_variable < other.m_variable ||
                 (m_variable == other.m_variable &&
                   m_symbol < other.m_symbol))));
    }
};

//--- start generated class pbes_equation ---//
/// \\brief list of pbes_equations
typedef atermpp::term_list<pbes_equation> pbes_equation_list;

/// \\brief vector of pbes_equations
typedef std::vector<pbes_equation>    pbes_equation_vector;

// prototype declaration
std::string pp(const pbes_equation& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const pbes_equation& x)
{
  return out << pbes_system::pp(x);
}

/// \\brief swap overload
inline void swap(pbes_equation& t1, pbes_equation& t2)
{
  t1.swap(t2);
}
//--- end generated class pbes_equation ---//

inline bool
operator==(const pbes_equation& x, const pbes_equation& y)
{
  return x.symbol() == y.symbol() &&
         x.variable() == y.variable() &&
         x.formula() == y.formula();
}

inline bool
operator!=(const pbes_equation& x, const pbes_equation& y)
{
  return !(x == y);
}

/// \brief Conversion to atermaPpl.
/// \return The specification converted to aterm format.
inline
atermpp::aterm pbes_equation_to_aterm(const pbes_equation& eqn)
{
  return atermpp::aterm(core::detail::function_symbol_PBEqn(), eqn.symbol(), eqn.variable(), eqn.formula());
}

// template function overloads
std::string pp(const pbes_equation_vector& x);
void normalize_sorts(pbes_equation_vector& x, const data::sort_specification& sortspec);
std::set<data::variable> find_free_variables(const pbes_system::pbes_equation& x);

} // namespace pbes_system

} // namespace mcrl2

#endif // MCRL2_PBES_PBES_EQUATION_H
