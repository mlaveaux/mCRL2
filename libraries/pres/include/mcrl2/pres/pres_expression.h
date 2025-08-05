// Author(s): Jan Friso Groote. Based on pbes_expression.h by Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pres/pres_expression.h
/// \brief The class pres_expression.

#ifndef MCRL2_PRES_PRES_EXPRESSION_H
#define MCRL2_PRES_PRES_EXPRESSION_H

#include "mcrl2/data/expression_traits.h"
#include "mcrl2/data/optimized_boolean_operators.h"
#include "mcrl2/data/real_utilities.h"
#include "mcrl2/data/cardinality.h"
#include "mcrl2/pbes/propositional_variable.h"

namespace mcrl2::pres_system
{

/// \brief The propositional variable is taken from a pbes_system.
using propositional_variable = pbes_system::propositional_variable;

//--- start generated classes ---//
/// \\brief A pres expression
class pres_expression: public atermpp::aterm
{
  public:
    /// \\brief Default constructor X3.
    pres_expression()
      : atermpp::aterm(core::detail::default_values::PRExpr)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit pres_expression(const atermpp::aterm& term)
      : atermpp::aterm(term)
    {
      assert(core::detail::check_rule_PRExpr(*this));
    }

    /// \\brief Constructor Z6.
    pres_expression(const data::data_expression& x)
      : atermpp::aterm(x)
    {}

    /// \\brief Constructor Z6.
    pres_expression(const data::variable& x)
      : atermpp::aterm(x)
    {}

    /// \\brief Constructor Z6.
    pres_expression(const data::untyped_data_parameter& x)
      : atermpp::aterm(x)
    {}

    /// Move semantics
    pres_expression(const pres_expression&) noexcept = default;
    pres_expression(pres_expression&&) noexcept = default;
    pres_expression& operator=(const pres_expression&) noexcept = default;
    pres_expression& operator=(pres_expression&&) noexcept = default;
};

/// \\brief list of pres_expressions
using pres_expression_list = atermpp::term_list<pres_expression>;

/// \\brief vector of pres_expressions
using pres_expression_vector = std::vector<pres_expression>;

// prototypes
inline bool is_propositional_variable_instantiation(const atermpp::aterm& x);
inline bool is_minus(const atermpp::aterm& x);
inline bool is_and(const atermpp::aterm& x);
inline bool is_or(const atermpp::aterm& x);
inline bool is_imp(const atermpp::aterm& x);
inline bool is_plus(const atermpp::aterm& x);
inline bool is_const_multiply(const atermpp::aterm& x);
inline bool is_const_multiply_alt(const atermpp::aterm& x);
inline bool is_infimum(const atermpp::aterm& x);
inline bool is_supremum(const atermpp::aterm& x);
inline bool is_sum(const atermpp::aterm& x);
inline bool is_eqinf(const atermpp::aterm& x);
inline bool is_eqninf(const atermpp::aterm& x);
inline bool is_condsm(const atermpp::aterm& x);
inline bool is_condeq(const atermpp::aterm& x);

/// \\brief Test for a pres_expression expression
/// \\param x A term
/// \\return True if \\a x is a pres_expression expression
inline
bool is_pres_expression(const atermpp::aterm& x)
{
  return data::is_data_expression(x) ||
         data::is_variable(x) ||
         data::is_untyped_data_parameter(x) ||
         pres_system::is_propositional_variable_instantiation(x) ||
         pres_system::is_minus(x) ||
         pres_system::is_and(x) ||
         pres_system::is_or(x) ||
         pres_system::is_imp(x) ||
         pres_system::is_plus(x) ||
         pres_system::is_const_multiply(x) ||
         pres_system::is_const_multiply_alt(x) ||
         pres_system::is_infimum(x) ||
         pres_system::is_supremum(x) ||
         pres_system::is_sum(x) ||
         pres_system::is_eqinf(x) ||
         pres_system::is_eqninf(x) ||
         pres_system::is_condsm(x) ||
         pres_system::is_condeq(x);
}

// prototype declaration
std::string pp(const pres_expression& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const pres_expression& x)
{
  return out << pres_system::pp(x);
}

/// \\brief swap overload
inline void swap(pres_expression& t1, pres_expression& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief A propositional variable instantiation
class propositional_variable_instantiation: public pres_expression
{
  public:


    /// Move semantics
    propositional_variable_instantiation(const propositional_variable_instantiation&) noexcept = default;
    propositional_variable_instantiation(propositional_variable_instantiation&&) noexcept = default;
    propositional_variable_instantiation& operator=(const propositional_variable_instantiation&) noexcept = default;
    propositional_variable_instantiation& operator=(propositional_variable_instantiation&&) noexcept = default;

    const core::identifier_string& name() const
    {
      return atermpp::down_cast<core::identifier_string>((*this)[0]);
    }

    const data::data_expression_list& parameters() const
    {
      return atermpp::down_cast<data::data_expression_list>((*this)[1]);
    }
//--- start user section propositional_variable_instantiation ---//
    /// \brief Default constructor.
    propositional_variable_instantiation()
      : pres_expression(core::detail::default_values::PropVarInst)
    {}

    /// \brief Constructor.
    /// \param term A term
    explicit propositional_variable_instantiation(const atermpp::aterm& term)
      : pres_expression(term)
    {
      assert(core::detail::check_term_PropVarInst(*this));
    }

    /// \brief Constructor.
    propositional_variable_instantiation(const core::identifier_string& name, const data::data_expression_list& parameters)
    {
      atermpp::make_term_appl(*this,core::detail::function_symbol_PropVarInst(), name, parameters);
    }

    /// \brief Constructor.
    propositional_variable_instantiation(const std::string& name, const data::data_expression_list& parameters)
    {
      atermpp::make_term_appl(*this,core::detail::function_symbol_PropVarInst(), core::identifier_string(name), parameters);
    }

//--- end user section propositional_variable_instantiation ---//
};

/// \\brief Make_propositional_variable_instantiation constructs a new term into a given address.
/// \\ \param t The reference into which the new propositional_variable_instantiation is constructed. 
template <class... ARGUMENTS>
inline void make_propositional_variable_instantiation(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_PropVarInst(), args...);
}

/// \\brief list of propositional_variable_instantiations
using propositional_variable_instantiation_list = atermpp::term_list<propositional_variable_instantiation>;

/// \\brief vector of propositional_variable_instantiations
using propositional_variable_instantiation_vector = std::vector<propositional_variable_instantiation>;

/// \\brief Test for a propositional_variable_instantiation expression
/// \\param x A term
/// \\return True if \\a x is a propositional_variable_instantiation expression
inline
bool is_propositional_variable_instantiation(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::PropVarInst;
}

// prototype declaration
std::string pp(const propositional_variable_instantiation& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const propositional_variable_instantiation& x)
{
  return out << pres_system::pp(x);
}

/// \\brief swap overload
inline void swap(propositional_variable_instantiation& t1, propositional_variable_instantiation& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The not operator for pres expressions
class minus: public pres_expression
{
  public:
    /// \\brief Default constructor X3.
    minus()
      : pres_expression(core::detail::default_values::PRESMinus)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit minus(const atermpp::aterm& term)
      : pres_expression(term)
    {
      assert(core::detail::check_term_PRESMinus(*this));
    }

    /// \\brief Constructor Z14.
    explicit minus(const pres_expression& operand)
      : pres_expression(atermpp::aterm(core::detail::function_symbol_PRESMinus(), operand))
    {}

    /// Move semantics
    minus(const minus&) noexcept = default;
    minus(minus&&) noexcept = default;
    minus& operator=(const minus&) noexcept = default;
    minus& operator=(minus&&) noexcept = default;

    const pres_expression& operand() const
    {
      return atermpp::down_cast<pres_expression>((*this)[0]);
    }
};

/// \\brief Make_minus constructs a new term into a given address.
/// \\ \param t The reference into which the new minus is constructed. 
template <class... ARGUMENTS>
inline void make_minus(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_PRESMinus(), args...);
}

/// \\brief Test for a minus expression
/// \\param x A term
/// \\return True if \\a x is a minus expression
inline
bool is_minus(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::PRESMinus;
}

// prototype declaration
std::string pp(const minus& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const minus& x)
{
  return out << pres_system::pp(x);
}

/// \\brief swap overload
inline void swap(minus& t1, minus& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The and operator for pres expressions
class and_: public pres_expression
{
  public:
    /// \\brief Default constructor X3.
    and_()
      : pres_expression(core::detail::default_values::PRESAnd)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit and_(const atermpp::aterm& term)
      : pres_expression(term)
    {
      assert(core::detail::check_term_PRESAnd(*this));
    }

    /// \\brief Constructor Z14.
    and_(const pres_expression& left, const pres_expression& right)
      : pres_expression(atermpp::aterm(core::detail::function_symbol_PRESAnd(), left, right))
    {}

    /// Move semantics
    and_(const and_&) noexcept = default;
    and_(and_&&) noexcept = default;
    and_& operator=(const and_&) noexcept = default;
    and_& operator=(and_&&) noexcept = default;

    const pres_expression& left() const
    {
      return atermpp::down_cast<pres_expression>((*this)[0]);
    }

    const pres_expression& right() const
    {
      return atermpp::down_cast<pres_expression>((*this)[1]);
    }
};

/// \\brief Make_and_ constructs a new term into a given address.
/// \\ \param t The reference into which the new and_ is constructed. 
template <class... ARGUMENTS>
inline void make_and_(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_PRESAnd(), args...);
}

/// \\brief Test for a and expression
/// \\param x A term
/// \\return True if \\a x is a and expression
inline
bool is_and(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::PRESAnd;
}

// prototype declaration
std::string pp(const and_& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const and_& x)
{
  return out << pres_system::pp(x);
}

/// \\brief swap overload
inline void swap(and_& t1, and_& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The or operator for pres expressions
class or_: public pres_expression
{
  public:
    /// \\brief Default constructor X3.
    or_()
      : pres_expression(core::detail::default_values::PRESOr)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit or_(const atermpp::aterm& term)
      : pres_expression(term)
    {
      assert(core::detail::check_term_PRESOr(*this));
    }

    /// \\brief Constructor Z14.
    or_(const pres_expression& left, const pres_expression& right)
      : pres_expression(atermpp::aterm(core::detail::function_symbol_PRESOr(), left, right))
    {}

    /// Move semantics
    or_(const or_&) noexcept = default;
    or_(or_&&) noexcept = default;
    or_& operator=(const or_&) noexcept = default;
    or_& operator=(or_&&) noexcept = default;

    const pres_expression& left() const
    {
      return atermpp::down_cast<pres_expression>((*this)[0]);
    }

    const pres_expression& right() const
    {
      return atermpp::down_cast<pres_expression>((*this)[1]);
    }
};

/// \\brief Make_or_ constructs a new term into a given address.
/// \\ \param t The reference into which the new or_ is constructed. 
template <class... ARGUMENTS>
inline void make_or_(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_PRESOr(), args...);
}

/// \\brief Test for a or expression
/// \\param x A term
/// \\return True if \\a x is a or expression
inline
bool is_or(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::PRESOr;
}

// prototype declaration
std::string pp(const or_& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const or_& x)
{
  return out << pres_system::pp(x);
}

/// \\brief swap overload
inline void swap(or_& t1, or_& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The implication operator for pres expressions
class imp: public pres_expression
{
  public:
    /// \\brief Default constructor X3.
    imp()
      : pres_expression(core::detail::default_values::PRESImp)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit imp(const atermpp::aterm& term)
      : pres_expression(term)
    {
      assert(core::detail::check_term_PRESImp(*this));
    }

    /// \\brief Constructor Z14.
    imp(const pres_expression& left, const pres_expression& right)
      : pres_expression(atermpp::aterm(core::detail::function_symbol_PRESImp(), left, right))
    {}

    /// Move semantics
    imp(const imp&) noexcept = default;
    imp(imp&&) noexcept = default;
    imp& operator=(const imp&) noexcept = default;
    imp& operator=(imp&&) noexcept = default;

    const pres_expression& left() const
    {
      return atermpp::down_cast<pres_expression>((*this)[0]);
    }

    const pres_expression& right() const
    {
      return atermpp::down_cast<pres_expression>((*this)[1]);
    }
};

/// \\brief Make_imp constructs a new term into a given address.
/// \\ \param t The reference into which the new imp is constructed. 
template <class... ARGUMENTS>
inline void make_imp(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_PRESImp(), args...);
}

/// \\brief Test for a imp expression
/// \\param x A term
/// \\return True if \\a x is a imp expression
inline
bool is_imp(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::PRESImp;
}

// prototype declaration
std::string pp(const imp& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const imp& x)
{
  return out << pres_system::pp(x);
}

/// \\brief swap overload
inline void swap(imp& t1, imp& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The addition operator for pres expressions
class plus: public pres_expression
{
  public:
    /// \\brief Default constructor X3.
    plus()
      : pres_expression(core::detail::default_values::PRESPlus)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit plus(const atermpp::aterm& term)
      : pres_expression(term)
    {
      assert(core::detail::check_term_PRESPlus(*this));
    }

    /// \\brief Constructor Z14.
    plus(const pres_expression& left, const pres_expression& right)
      : pres_expression(atermpp::aterm(core::detail::function_symbol_PRESPlus(), left, right))
    {}

    /// Move semantics
    plus(const plus&) noexcept = default;
    plus(plus&&) noexcept = default;
    plus& operator=(const plus&) noexcept = default;
    plus& operator=(plus&&) noexcept = default;

    const pres_expression& left() const
    {
      return atermpp::down_cast<pres_expression>((*this)[0]);
    }

    const pres_expression& right() const
    {
      return atermpp::down_cast<pres_expression>((*this)[1]);
    }
};

/// \\brief Make_plus constructs a new term into a given address.
/// \\ \param t The reference into which the new plus is constructed. 
template <class... ARGUMENTS>
inline void make_plus(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_PRESPlus(), args...);
}

/// \\brief Test for a plus expression
/// \\param x A term
/// \\return True if \\a x is a plus expression
inline
bool is_plus(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::PRESPlus;
}

// prototype declaration
std::string pp(const plus& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const plus& x)
{
  return out << pres_system::pp(x);
}

/// \\brief swap overload
inline void swap(plus& t1, plus& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The multiplication with a positive constant with the constant at the left.
class const_multiply: public pres_expression
{
  public:
    /// \\brief Default constructor X3.
    const_multiply()
      : pres_expression(core::detail::default_values::PRESConstantMultiply)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit const_multiply(const atermpp::aterm& term)
      : pres_expression(term)
    {
      assert(core::detail::check_term_PRESConstantMultiply(*this));
    }

    /// \\brief Constructor Z14.
    const_multiply(const data::data_expression& left, const pres_expression& right)
      : pres_expression(atermpp::aterm(core::detail::function_symbol_PRESConstantMultiply(), left, right))
    {}

    /// Move semantics
    const_multiply(const const_multiply&) noexcept = default;
    const_multiply(const_multiply&&) noexcept = default;
    const_multiply& operator=(const const_multiply&) noexcept = default;
    const_multiply& operator=(const_multiply&&) noexcept = default;

    const data::data_expression& left() const
    {
      return atermpp::down_cast<data::data_expression>((*this)[0]);
    }

    const pres_expression& right() const
    {
      return atermpp::down_cast<pres_expression>((*this)[1]);
    }
};

/// \\brief Make_const_multiply constructs a new term into a given address.
/// \\ \param t The reference into which the new const_multiply is constructed. 
template <class... ARGUMENTS>
inline void make_const_multiply(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_PRESConstantMultiply(), args...);
}

/// \\brief Test for a const_multiply expression
/// \\param x A term
/// \\return True if \\a x is a const_multiply expression
inline
bool is_const_multiply(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::PRESConstantMultiply;
}

// prototype declaration
std::string pp(const const_multiply& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const const_multiply& x)
{
  return out << pres_system::pp(x);
}

/// \\brief swap overload
inline void swap(const_multiply& t1, const_multiply& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The multiplication with a positive constant with the constant at the right.
class const_multiply_alt: public pres_expression
{
  public:
    /// \\brief Default constructor X3.
    const_multiply_alt()
      : pres_expression(core::detail::default_values::PRESConstantMultiplyAlt)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit const_multiply_alt(const atermpp::aterm& term)
      : pres_expression(term)
    {
      assert(core::detail::check_term_PRESConstantMultiplyAlt(*this));
    }

    /// \\brief Constructor Z14.
    const_multiply_alt(const pres_expression& left, const data::data_expression& right)
      : pres_expression(atermpp::aterm(core::detail::function_symbol_PRESConstantMultiplyAlt(), left, right))
    {}

    /// Move semantics
    const_multiply_alt(const const_multiply_alt&) noexcept = default;
    const_multiply_alt(const_multiply_alt&&) noexcept = default;
    const_multiply_alt& operator=(const const_multiply_alt&) noexcept = default;
    const_multiply_alt& operator=(const_multiply_alt&&) noexcept = default;

    const pres_expression& left() const
    {
      return atermpp::down_cast<pres_expression>((*this)[0]);
    }

    const data::data_expression& right() const
    {
      return atermpp::down_cast<data::data_expression>((*this)[1]);
    }
};

/// \\brief Make_const_multiply_alt constructs a new term into a given address.
/// \\ \param t The reference into which the new const_multiply_alt is constructed. 
template <class... ARGUMENTS>
inline void make_const_multiply_alt(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_PRESConstantMultiplyAlt(), args...);
}

/// \\brief Test for a const_multiply_alt expression
/// \\param x A term
/// \\return True if \\a x is a const_multiply_alt expression
inline
bool is_const_multiply_alt(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::PRESConstantMultiplyAlt;
}

// prototype declaration
std::string pp(const const_multiply_alt& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const const_multiply_alt& x)
{
  return out << pres_system::pp(x);
}

/// \\brief swap overload
inline void swap(const_multiply_alt& t1, const_multiply_alt& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The infimum over a data type for pres expressions
class infimum: public pres_expression
{
  public:
    /// \\brief Default constructor X3.
    infimum()
      : pres_expression(core::detail::default_values::PRESInfimum)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit infimum(const atermpp::aterm& term)
      : pres_expression(term)
    {
      assert(core::detail::check_term_PRESInfimum(*this));
    }

    /// \\brief Constructor Z14.
    infimum(const data::variable_list& variables, const pres_expression& body)
      : pres_expression(atermpp::aterm(core::detail::function_symbol_PRESInfimum(), variables, body))
    {}

    /// Move semantics
    infimum(const infimum&) noexcept = default;
    infimum(infimum&&) noexcept = default;
    infimum& operator=(const infimum&) noexcept = default;
    infimum& operator=(infimum&&) noexcept = default;

    const data::variable_list& variables() const
    {
      return atermpp::down_cast<data::variable_list>((*this)[0]);
    }

    const pres_expression& body() const
    {
      return atermpp::down_cast<pres_expression>((*this)[1]);
    }
};

/// \\brief Make_infimum constructs a new term into a given address.
/// \\ \param t The reference into which the new infimum is constructed. 
template <class... ARGUMENTS>
inline void make_infimum(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_PRESInfimum(), args...);
}

/// \\brief Test for a infimum expression
/// \\param x A term
/// \\return True if \\a x is a infimum expression
inline
bool is_infimum(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::PRESInfimum;
}

// prototype declaration
std::string pp(const infimum& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const infimum& x)
{
  return out << pres_system::pp(x);
}

/// \\brief swap overload
inline void swap(infimum& t1, infimum& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The supremeum over a data type for pres expressions
class supremum: public pres_expression
{
  public:
    /// \\brief Default constructor X3.
    supremum()
      : pres_expression(core::detail::default_values::PRESSupremum)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit supremum(const atermpp::aterm& term)
      : pres_expression(term)
    {
      assert(core::detail::check_term_PRESSupremum(*this));
    }

    /// \\brief Constructor Z14.
    supremum(const data::variable_list& variables, const pres_expression& body)
      : pres_expression(atermpp::aterm(core::detail::function_symbol_PRESSupremum(), variables, body))
    {}

    /// Move semantics
    supremum(const supremum&) noexcept = default;
    supremum(supremum&&) noexcept = default;
    supremum& operator=(const supremum&) noexcept = default;
    supremum& operator=(supremum&&) noexcept = default;

    const data::variable_list& variables() const
    {
      return atermpp::down_cast<data::variable_list>((*this)[0]);
    }

    const pres_expression& body() const
    {
      return atermpp::down_cast<pres_expression>((*this)[1]);
    }
};

/// \\brief Make_supremum constructs a new term into a given address.
/// \\ \param t The reference into which the new supremum is constructed. 
template <class... ARGUMENTS>
inline void make_supremum(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_PRESSupremum(), args...);
}

/// \\brief Test for a supremum expression
/// \\param x A term
/// \\return True if \\a x is a supremum expression
inline
bool is_supremum(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::PRESSupremum;
}

// prototype declaration
std::string pp(const supremum& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const supremum& x)
{
  return out << pres_system::pp(x);
}

/// \\brief swap overload
inline void swap(supremum& t1, supremum& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The generic sum operator for pres expressions
class sum: public pres_expression
{
  public:
    /// \\brief Default constructor X3.
    sum()
      : pres_expression(core::detail::default_values::PRESSum)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit sum(const atermpp::aterm& term)
      : pres_expression(term)
    {
      assert(core::detail::check_term_PRESSum(*this));
    }

    /// \\brief Constructor Z14.
    sum(const data::variable_list& variables, const pres_expression& body)
      : pres_expression(atermpp::aterm(core::detail::function_symbol_PRESSum(), variables, body))
    {}

    /// Move semantics
    sum(const sum&) noexcept = default;
    sum(sum&&) noexcept = default;
    sum& operator=(const sum&) noexcept = default;
    sum& operator=(sum&&) noexcept = default;

    const data::variable_list& variables() const
    {
      return atermpp::down_cast<data::variable_list>((*this)[0]);
    }

    const pres_expression& body() const
    {
      return atermpp::down_cast<pres_expression>((*this)[1]);
    }
};

/// \\brief Make_sum constructs a new term into a given address.
/// \\ \param t The reference into which the new sum is constructed. 
template <class... ARGUMENTS>
inline void make_sum(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_PRESSum(), args...);
}

/// \\brief Test for a sum expression
/// \\param x A term
/// \\return True if \\a x is a sum expression
inline
bool is_sum(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::PRESSum;
}

// prototype declaration
std::string pp(const sum& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const sum& x)
{
  return out << pres_system::pp(x);
}

/// \\brief swap overload
inline void swap(sum& t1, sum& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The indicator whether the argument is infinite
class eqinf: public pres_expression
{
  public:
    /// \\brief Default constructor X3.
    eqinf()
      : pres_expression(core::detail::default_values::PRESEqInf)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit eqinf(const atermpp::aterm& term)
      : pres_expression(term)
    {
      assert(core::detail::check_term_PRESEqInf(*this));
    }

    /// \\brief Constructor Z14.
    explicit eqinf(const pres_expression& operand)
      : pres_expression(atermpp::aterm(core::detail::function_symbol_PRESEqInf(), operand))
    {}

    /// Move semantics
    eqinf(const eqinf&) noexcept = default;
    eqinf(eqinf&&) noexcept = default;
    eqinf& operator=(const eqinf&) noexcept = default;
    eqinf& operator=(eqinf&&) noexcept = default;

    const pres_expression& operand() const
    {
      return atermpp::down_cast<pres_expression>((*this)[0]);
    }
};

/// \\brief Make_eqinf constructs a new term into a given address.
/// \\ \param t The reference into which the new eqinf is constructed. 
template <class... ARGUMENTS>
inline void make_eqinf(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_PRESEqInf(), args...);
}

/// \\brief Test for a eqinf expression
/// \\param x A term
/// \\return True if \\a x is a eqinf expression
inline
bool is_eqinf(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::PRESEqInf;
}

// prototype declaration
std::string pp(const eqinf& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const eqinf& x)
{
  return out << pres_system::pp(x);
}

/// \\brief swap overload
inline void swap(eqinf& t1, eqinf& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The indicator whether the argument is -infinite
class eqninf: public pres_expression
{
  public:
    /// \\brief Default constructor X3.
    eqninf()
      : pres_expression(core::detail::default_values::PRESEqNInf)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit eqninf(const atermpp::aterm& term)
      : pres_expression(term)
    {
      assert(core::detail::check_term_PRESEqNInf(*this));
    }

    /// \\brief Constructor Z14.
    explicit eqninf(const pres_expression& operand)
      : pres_expression(atermpp::aterm(core::detail::function_symbol_PRESEqNInf(), operand))
    {}

    /// Move semantics
    eqninf(const eqninf&) noexcept = default;
    eqninf(eqninf&&) noexcept = default;
    eqninf& operator=(const eqninf&) noexcept = default;
    eqninf& operator=(eqninf&&) noexcept = default;

    const pres_expression& operand() const
    {
      return atermpp::down_cast<pres_expression>((*this)[0]);
    }
};

/// \\brief Make_eqninf constructs a new term into a given address.
/// \\ \param t The reference into which the new eqninf is constructed. 
template <class... ARGUMENTS>
inline void make_eqninf(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_PRESEqNInf(), args...);
}

/// \\brief Test for a eqninf expression
/// \\param x A term
/// \\return True if \\a x is a eqninf expression
inline
bool is_eqninf(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::PRESEqNInf;
}

// prototype declaration
std::string pp(const eqninf& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const eqninf& x)
{
  return out << pres_system::pp(x);
}

/// \\brief swap overload
inline void swap(eqninf& t1, eqninf& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief Conditional operator with condition smaller than 0
class condsm: public pres_expression
{
  public:
    /// \\brief Default constructor X3.
    condsm()
      : pres_expression(core::detail::default_values::PRESCondSm)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit condsm(const atermpp::aterm& term)
      : pres_expression(term)
    {
      assert(core::detail::check_term_PRESCondSm(*this));
    }

    /// \\brief Constructor Z14.
    condsm(const pres_expression& arg1, const pres_expression& arg2, const pres_expression& arg3)
      : pres_expression(atermpp::aterm(core::detail::function_symbol_PRESCondSm(), arg1, arg2, arg3))
    {}

    /// Move semantics
    condsm(const condsm&) noexcept = default;
    condsm(condsm&&) noexcept = default;
    condsm& operator=(const condsm&) noexcept = default;
    condsm& operator=(condsm&&) noexcept = default;

    const pres_expression& arg1() const
    {
      return atermpp::down_cast<pres_expression>((*this)[0]);
    }

    const pres_expression& arg2() const
    {
      return atermpp::down_cast<pres_expression>((*this)[1]);
    }

    const pres_expression& arg3() const
    {
      return atermpp::down_cast<pres_expression>((*this)[2]);
    }
};

/// \\brief Make_condsm constructs a new term into a given address.
/// \\ \param t The reference into which the new condsm is constructed. 
template <class... ARGUMENTS>
inline void make_condsm(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_PRESCondSm(), args...);
}

/// \\brief Test for a condsm expression
/// \\param x A term
/// \\return True if \\a x is a condsm expression
inline
bool is_condsm(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::PRESCondSm;
}

// prototype declaration
std::string pp(const condsm& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const condsm& x)
{
  return out << pres_system::pp(x);
}

/// \\brief swap overload
inline void swap(condsm& t1, condsm& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief Conditional operator with condition smaller than or equal to 0
class condeq: public pres_expression
{
  public:
    /// \\brief Default constructor X3.
    condeq()
      : pres_expression(core::detail::default_values::PRESCondEq)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit condeq(const atermpp::aterm& term)
      : pres_expression(term)
    {
      assert(core::detail::check_term_PRESCondEq(*this));
    }

    /// \\brief Constructor Z14.
    condeq(const pres_expression& arg1, const pres_expression& arg2, const pres_expression& arg3)
      : pres_expression(atermpp::aterm(core::detail::function_symbol_PRESCondEq(), arg1, arg2, arg3))
    {}

    /// Move semantics
    condeq(const condeq&) noexcept = default;
    condeq(condeq&&) noexcept = default;
    condeq& operator=(const condeq&) noexcept = default;
    condeq& operator=(condeq&&) noexcept = default;

    const pres_expression& arg1() const
    {
      return atermpp::down_cast<pres_expression>((*this)[0]);
    }

    const pres_expression& arg2() const
    {
      return atermpp::down_cast<pres_expression>((*this)[1]);
    }

    const pres_expression& arg3() const
    {
      return atermpp::down_cast<pres_expression>((*this)[2]);
    }
};

/// \\brief Make_condeq constructs a new term into a given address.
/// \\ \param t The reference into which the new condeq is constructed. 
template <class... ARGUMENTS>
inline void make_condeq(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_PRESCondEq(), args...);
}

/// \\brief Test for a condeq expression
/// \\param x A term
/// \\return True if \\a x is a condeq expression
inline
bool is_condeq(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::PRESCondEq;
}

// prototype declaration
std::string pp(const condeq& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const condeq& x)
{
  return out << pres_system::pp(x);
}

/// \\brief swap overload
inline void swap(condeq& t1, condeq& t2) noexcept
{
  t1.swap(t2);
}
//--- end generated classes ---//

// template function overloads
std::string pp(const pres_expression_list& x);
std::string pp(const pres_expression_vector& x);
std::string pp(const propositional_variable_instantiation_list& x);
std::string pp(const propositional_variable_instantiation_vector& x);
std::set<pres_system::propositional_variable_instantiation> find_propositional_variable_instantiations(const pres_system::pres_expression& x);
std::set<core::identifier_string> find_identifiers(const pres_system::pres_expression& x);
std::set<data::variable> find_free_variables(const pres_system::pres_expression& x);
bool search_variable(const pres_system::pres_expression& x, const data::variable& v);
pres_system::pres_expression normalize_sorts(const pres_system::pres_expression& x, const data::sort_specification& sortspec);
pres_system::pres_expression translate_user_notation(const pres_system::pres_expression& x);

/// \return Returns the value true
inline
const pres_expression& true_()
{
  /* The dynamic cast is required, to prevent copying the data term true to
     a local term on the stack. */
  return reinterpret_cast<const pres_expression&>(data::sort_bool::true_());
}

/// \return Returns the value false
inline
const pres_expression& false_()
{
  /* The dynamic cast is required, to prevent copying the data term false to
     a local term on the stack. */
  return reinterpret_cast<const pres_expression&>(data::sort_bool::false_());
}

/// \brief Test for the value true
/// \param t A PRES expression
/// \return True if it is the value \p true
inline bool is_true(const pres_expression& t)
{
  return data::sort_bool::is_true_function_symbol(t);
}

/// \brief Test for the value false
/// \param t A PRES expression
/// \return True if it is the value \p false
inline bool is_false(const pres_expression& t)
{
  return data::sort_bool::is_false_function_symbol(t);
}

/// \brief Returns true if the term t is a minus expression
/// \param t A PRES expression
/// \return True if the term t is a minus expression
inline bool is_pres_minus(const pres_expression& t)
{
  return pres_system::is_minus(t);
}

/// \brief Returns true if the term t is an and expression
/// \param t A PRES expression
/// \return True if the term t is an and expression
inline bool is_pres_and(const pres_expression& t)
{
  return pres_system::is_and(t);
}

/// \brief Returns true if the term t is an or expression
/// \param t A PRES expression
/// \return True if the term t is an or expression
inline bool is_pres_or(const pres_expression& t)
{
  return pres_system::is_or(t);
}

/// \brief Returns true if the term t is an imp expression
/// \param t A PRES expression
/// \return True if the term t is an imp expression
inline bool is_pres_imp(const pres_expression& t)
{
  return pres_system::is_imp(t);
}

/// \brief Returns true if the term t is a generalized minus expression
/// \param t A PRES expression
/// \return True if the term t is a generalized minus expression
inline bool is_pres_infimum(const pres_expression& t)
{
  return pres_system::is_infimum(t);
}

/// \brief Returns true if the term t is a generalized maximum expression
/// \param t A PRES expression
/// \return True if the term t is a generalized maximum expression
inline bool is_pres_supremum(const pres_expression& t)
{
  return pres_system::is_supremum(t);
}

/* /// \brief Test for a conjunction
/// \param t A PRES expression or a data expression
/// \return True if it is a conjunction
inline bool is_universal_not(const pres_expression& t)
{
  return is_pres_not(t) || data::sort_bool::is_not_application(t);
} */

/// \brief Test for a conjunction
/// \param t A PRES expression or a data expression
/// \return True if it is a conjunction
inline bool is_universal_and(const pres_expression& t)
{
  return is_pres_and(t) || data::sort_bool::is_and_application(t);
}

/// \brief Test for a disjunction
/// \param t A PRES expression or a data expression
/// \return True if it is a disjunction
inline bool is_universal_or(const pres_expression& t)
{
  return is_pres_or(t) || data::sort_bool::is_or_application(t);
}

/// \brief Returns true if the term t is a data expression
/// \param t A PRES expression
/// \return True if the term t is a data expression
inline bool is_data(const pres_expression& t)
{
  return data::is_data_expression(t);
}

/// \brief The namespace for accessor functions on pres expressions.
namespace accessors
{

/// \brief Returns the pres expression argument of expressions of type not, exists and forall.
/// \param t A PRES expression
/// \return The pres expression argument of expressions of type not, exists and forall.
inline
const pres_expression& arg(const pres_expression& t)
{
  if (is_pres_minus(t))
  {
    return atermpp::down_cast<const pres_expression>(t[0]);
  }
  else
  {
    assert(is_infimum(t) || is_supremum(t) || is_sum(t));
    return atermpp::down_cast<const pres_expression>(t[1]);
  }
}

/// \brief Returns the pres expression argument of expressions of type not, exists and forall.
/// \param t A PRES expression or a data expression
/// \return The pres expression argument of expressions of type not, exists and forall.
inline
pres_expression data_arg(const pres_expression& t)
{
  if (data::is_data_expression(t))
  {
    assert(data::is_application(t));
    const auto& a = atermpp::down_cast<const data::application>(t);
    return *(a.begin());
  }
  else
  {
    return arg(t);
  }
}

/// \brief Returns the left hand side of an expression of type and, or or imp.
/// \param t A PRES expression
/// \return The left hand side of an expression of type and, or or imp.
inline
const pres_expression& left(const pres_expression& t)
{
  assert(is_and(t) || is_or(t) || is_imp(t));
  return atermpp::down_cast<pres_expression>(t[0]);
}

/// \brief Returns the left hand side of an expression of type and, or or imp.
/// \param x A PRES expression or a data expression
/// \return The left hand side of an expression of type and, or or imp.
inline
pres_expression data_left(const pres_expression& x)
{
  if (data::is_data_expression(x))
  {
    return data::binary_left(atermpp::down_cast<data::application>(x));
  }
  else
  {
    return left(x);
  }
}

/// \brief Returns the right hand side of an expression of type and, or or imp.
/// \param t A PRES expression
/// \return The right hand side of an expression of type and, or or imp.
inline
const pres_expression& right(const pres_expression& t)
{
  return atermpp::down_cast<pres_expression>(t[1]);
}

/// \brief Returns the left hand side of an expression of type and, or or imp.
/// \param x A PRES expression or a data expression
/// \return The left hand side of an expression of type and, or or imp.
inline
pres_expression data_right(const pres_expression& x)
{
  if (data::is_data_expression(x))
  {
    return data::binary_right(atermpp::down_cast<data::application>(x));
  }
  else
  {
    return right(x);
  }
}

/// \brief Returns the variables of a quantification expression
/// \param t A PRES expression
/// \return The variables of a quantification expression
inline
const data::variable_list& var(const pres_expression& t)
{
  assert(is_infimum(t) || is_supremum(t) || is_sum(t));
  return atermpp::down_cast<data::variable_list>(t[0]);
}

/// \brief Returns the name of a propositional variable expression
/// \param t A PRES expression
/// \return The name of a propositional variable expression
inline
const core::identifier_string& name(const pres_expression& t)
{
  assert(is_propositional_variable_instantiation(t));
  return atermpp::down_cast<core::identifier_string>(t[0]);
}

/// \brief Returns the parameters of a propositional variable instantiation.
/// \param t A PRES expression
/// \return The parameters of a propositional variable instantiation.
inline
const data::data_expression_list& param(const pres_expression& t)
{
  assert(is_propositional_variable_instantiation(t));
  return atermpp::down_cast<data::data_expression_list>(t[1]);
}
} // namespace accessors

/// \brief Make a generalized minimum. It checks for an empty variable list,
/// which is not allowed.
/// \param l A sequence of data variables
/// \param p A PRES expression
/// \return The value <tt>infimum l.p</tt>
inline
pres_expression make_infimum(const data::variable_list& l, const pres_expression& p)
{
  if (l.empty())
  {
    return p;
  }
  return pres_expression(atermpp::aterm(core::detail::function_symbol_PRESInfimum(), l, p));
}

/// \brief Make an generalized maximum. It checks for an empty variable list,
/// which is not allowed.
/// \param l A sequence of data variables
/// \param p A PRES expression
/// \return The value <tt>exists l.p</tt>
inline
pres_expression make_supremum(const data::variable_list& l, const pres_expression& p)
{
  if (l.empty())
  {
    return p;
  }
  return pres_expression(atermpp::aterm(core::detail::function_symbol_PRESSupremum(), l, p));
}

/// \brief Make a negation
/// \param p A PRES expression
/// \return The value <tt>!p</tt>
inline
void optimized_minus(pres_expression& result, const pres_expression& p)
{
  // Should be optimized. 
  // true false special cases.
  // data::sort_real::is_zero(p)
  make_minus(result, p);
}

/// \brief Make a conjunction
/// \param p A PRES expression
/// \param q A PRES expression
/// \return The value <tt>p && q</tt>
inline
void optimized_and(pres_expression& result, const pres_expression& p, const pres_expression& q)
{
  data::optimized_and(result, p, q);
}

/// \brief Make a disjunction
/// \param p A PRES expression
/// \param q A PRES expression
/// \return The value <tt>p || q</tt>
inline
void optimized_or(pres_expression& result, const pres_expression& p, const pres_expression& q)
{
  data::optimized_or(result, p, q);
}

/// \brief Make a disjunction
/// \param p A PRES expression
/// \param q A PRES expression
/// \return The value <tt>p || q</tt>
inline
void optimized_plus(pres_expression& result, const pres_expression& p, const pres_expression& q)
{
  if (is_true(p))
  {
    result=p;
  }
  else if (is_true(q))
  {
    result=q;
  }
  else if (is_false(p))
  {
    if (data::is_data_expression(q))
    {
      const data::data_expression& qa = atermpp::down_cast<data::data_expression>(q);
      if (qa.sort() == data::sort_real::real_())
      {
        result=p;
        return;
      }
    }
    make_plus(result, p, q);  
  }
  else if (is_false(q))
  {
    if (data::is_data_expression(p))
    {
      const data::data_expression& pa = atermpp::down_cast<data::data_expression>(p);
      if (pa.sort() == data::sort_real::real_())
      {
        result=q;
        return;
      }
    }
    make_plus(result, p, q);  
  }
  else if (data::sort_real::is_zero(p))
  {
    result=q;
  }
  else if (data::sort_real::is_zero(q))
  {
    result=p;
  }
  else
  {
    make_plus(result, p, q);
  }
}

/// \brief Make an implication
/// \param p A PRES expression
/// \param q A PRES expression
/// \return The value <tt>p => q</tt>
/* inline
void optimized_imp(pres_expression& result, const pres_expression& p, const pres_expression& q)
{
  data::optimized_imp(result, p, q);
} */

/// \brief Make an infimum quantification
/// If l is empty, p is returned.
/// \param l A sequence of data variables
/// \param p A PRES expression
/// \return The value <tt>inf l.p</tt>
inline
void optimized_infimum(pres_expression& result, const data::variable_list& l, const pres_expression& p)
{
  std::set<data::variable> free_variables = find_free_variables(p);
  data::variable_list new_l(l.begin(), 
                           l.end(), 
                           [](const data::variable& v){ return v; }, 
                           [&free_variables](const data::variable& v){ return free_variables.find(v)!=free_variables.end(); });
  if (new_l.empty())
  {
    result=p;
  }
  else
  { 
    make_infimum(result, new_l, p);
  }
}

/// \brief Make a supremum.
/// If l is empty, p is returned.
/// \param l A sequence of data variables
/// \param p A PRES expression
/// \return The value <tt>sup l.p</tt>
inline
void optimized_supremum(pres_expression& result, const data::variable_list& l, const pres_expression& p)
{
  std::set<data::variable> free_variables = find_free_variables(p);
  data::variable_list new_l(l.begin(), 
                           l.end(), 
                           [](const data::variable& v){ return v; }, 
                           [&free_variables](const data::variable& v){ return free_variables.find(v)!=free_variables.end(); });
  if (new_l.empty())
  {
    result=p;
  }
  else
  { 
    make_supremum(result, new_l, p);
  }
}

/// \brief Make an sum quantification.
/// If l is empty, p is returned.
/// \param l A sequence of data variables.
/// \param p A PRES expression.
/// \param data_specification A data specification to determine the cardinality of sorts.
/// \param rewr A rewriter to determine the cardinality of sorts.
/// \return The value <tt>sum l.p</tt>
inline
void optimized_sum(pres_expression& result, 
                   const data::variable_list& l, 
                   const pres_expression& p, 
                   const data::data_specification& data_specification, 
                   const data::rewriter& rewr)
{
  if (l.size()==0)
  {
    result=p;
    return;
  }
  std::set<data::variable> free_variables = find_free_variables(p);
  data::variable_list new_l;
  std::size_t factor=1;
  data::cardinality_calculator cardinality(data_specification, rewr);
  for(const data::variable& v: l)
  {
    if (free_variables.find(v)==free_variables.end()) // not found.
    {
      if (!is_true(p) && !is_false(p) && !data::sort_real::is_zero(p))
      {
        // Determine the cardinality. 
        std::size_t c = cardinality(v.sort());
        if (c==0) // This means the cardinality is infinite or cannot be determined.
        {
          new_l.push_front(v);
        }
        else
        {
          factor=factor*c;
        }
      }
    }
    else
    {
      new_l.push_front(v);
    }
  }
  if (factor!=1)
  {
    if (is_const_multiply(p))
    {
      const const_multiply& cm=atermpp::down_cast<const_multiply>(p);
      make_const_multiply(result, rewr(data::sort_real::times(cm.left(),data::sort_real::real_(factor))),cm.right());
    }
    else if (is_false(p) || is_true(p))
    {
      result=p;
    }
    else if (data::is_data_expression(p))
    {
      const data::data_expression& d = atermpp::down_cast<data::data_expression>(p);
      if (d.sort()==data::sort_bool::bool_())
      {
        result = d;  
      }
      else
      {
        assert(d.sort()==data::sort_real::real_());
        result=rewr(data::sort_real::times(d,data::sort_real::real_(factor)));
      }
    }
    else
    {
      make_const_multiply(result, data::sort_real::real_(factor),p);
    }
  }
  else 
  {
    result=p;
  }
  if (new_l.size()>0)
  {
    make_sum(result, new_l, result);
  }
  return;
} 

/// \brief Make an optimized condsm expression
/// \param p1 A pres expression
/// \param p2 A pres expression
/// \param p3 A pres expression
/// \return An optimized representation of condsm(p1, p2, p3)
inline
void optimized_condsm(pres_expression& result, const pres_expression& p1, const pres_expression& p2, const pres_expression& p3)
{
  if (p1==false_())
  {
    result = p2;
    return;
  }
  else if (p1==true_())
  {
    optimized_or(result, p2, p3);
    return;
  }
  make_condsm(result, p1, p2, p3);
  return;
}

/// \brief Make an optimized condsm expression
/// \param p1 A pres expression
/// \param p2 A pres expression
/// \param p3 A pres expression
/// \return An optimized representation of condsm(p1, p2, p3)
inline
void optimized_condeq(pres_expression& result, const pres_expression& p1, const pres_expression& p2, const pres_expression& p3)
{
  if (p1==false_())
  {
    // N.B. Here we use the fact that mCRL2 data types are never empty.
    optimized_and(result, p2, p3);
    return;
  }
  else if (p1==true_())
  {
    result = p3;
    return;
  }
  make_condeq(result, p1, p2, p3);
  return;
}

/// \brief Make an optimized eqinf expression
/// \param p A pres expression
/// \return An optimized representation of eqinf(p)
inline
void optimized_eqinf(pres_expression& result, const pres_expression& p)
{
  if (p==false_())
  {
    // N.B. Here we use the fact that mCRL2 data types are never empty.
    result = p;
    return;
  }
  else if (p==true_())
  {
    // N.B. Here we use the fact that mCRL2 data types are never empty.
    result = p;
    return;
  }
  else if (data::is_data_expression(p) && atermpp::down_cast<data::data_expression>(p).sort()==data::sort_real::real_())
  {
    result = false_();
    return;
  }
  make_eqinf(result, p);
}

/// \brief Make an optimized eqinf expression
/// \param p A pres expression
/// \return An optimized representation of eqinf(p)
inline
void optimized_eqninf(pres_expression& result, const pres_expression& p)
{
  if (p==false_())
  {
    // N.B. Here we use the fact that mCRL2 data types are never empty.
    result = p;
    return;
  }
  else if (p==true_())
  {
    // N.B. Here we use the fact that mCRL2 data types are never empty.
    result = p;
    return;
  }
  else if (data::is_data_expression(p) && atermpp::down_cast<data::data_expression>(p).sort()==data::sort_real::real_())
  {
    result = true_();
    return;
  }
  make_eqninf(result, p);
}

/// \brief Make an optimized const_multiply expression
/// \param d A real data expression, which should be positive. 
/// \param p A pres expression
/// \return An optimized representation of eqinf(p)
inline
void optimized_const_multiply(pres_expression& result, const data::data_expression& d, const pres_expression& p)
{
  if (data::sort_real::is_zero(d))
  {
    result = d;
    return;
  }
  if (data::sort_real::is_larger_zero(d) && 
      (p==false_() || p==true_() || is_eqinf(p) || is_eqninf(d)))
  {
    result = p;
    return;
  }
  make_const_multiply(result, d, p);
}

/// \brief Make an optimized const_multiply_alt expression
/// \param d A real data expression, which should be positive. 
/// \param p A pres expression
/// \return An optimized representation of eqinf(p)
inline
void optimized_const_multiply_alt(pres_expression& result, const data::data_expression& d, const pres_expression& p)
{
  if (data::sort_real::is_zero(d))
  { 
    result = d;
    return;
  }
  if (data::sort_real::is_larger_zero(d) && 
      (p==false_() || p==true_() || is_eqinf(p) || is_eqninf(d)))
  { 
    result = p;
    return;
  }
  make_const_multiply_alt(result, d, p);
}

inline
bool is_constant(const pres_expression& x)
{
  return find_free_variables(x).empty();
}

/* inline
const data::variable_list& quantifier_variables(const pres_expression& x)
{
  assert(is_exists(x) || is_forall(x));
  if (is_exists(x))
  {
    return atermpp::down_cast<exists>(x).variables();
  }
  else
  {
    return atermpp::down_cast<forall>(x).variables();
  }
} */

inline
data::variable_list free_variables(const pres_expression& x)
{
  std::set<data::variable> v = find_free_variables(x);
  return data::variable_list(v.begin(), v.end());
}

/// \\brief Make_propositional_variable constructs a new term into a given address.
/// \\ \param t The reference into which the new propositional_variable is constructed. 
template <class... ARGUMENTS>
inline void make_propositional_variable(atermpp::aterm& t, const ARGUMENTS&... args)
{
  mcrl2::pbes_system::make_propositional_variable(t, args...);
}

/// \\brief list of propositional_variables
using propositional_variable_list = atermpp::term_list<pbes_system::propositional_variable>;

/// \\brief vector of propositional_variables
using propositional_variable_vector = std::vector<pbes_system::propositional_variable>;

} // namespace mcrl2::pres_system

namespace mcrl2::core
{

/// \brief Contains type information for pres expressions.
template <>
struct term_traits<pres_system::pres_expression>
{
  /// \brief The term type
  using term_type = pres_system::pres_expression;

  /// \brief The data term type
  using data_term_type = data::data_expression;

  /// \brief The data term sequence type
  using data_term_sequence_type = data::data_expression_list;

  /// \brief The variable type
  using variable_type = data::variable;

  /// \brief The variable sequence type
  using variable_sequence_type = data::variable_list;

  /// \brief The propositional variable declaration type
  using propositional_variable_decl_type = pres_system::propositional_variable;

  /// \brief The propositional variable instantiation type
  using propositional_variable_type = pres_system::propositional_variable_instantiation;

  /// \brief The string type
  using string_type = core::identifier_string;

  /// \brief Make the value true
  /// \return The value \p true
  static inline
  term_type true_()
  {
    return pres_system::true_();
  }

  /// \brief Make the value false
  /// \return The value \p false
  static inline
  term_type false_()
  {
    return pres_system::false_();
  }

  /// \brief Make a negation
  /// \param p A term
  /// \return The value <tt>!p</tt>
  static inline
  term_type minus(const term_type& p)
  {
    return pres_system::minus(p);
  }

  /// \brief Make a negation
  /// \param result The value <tt>!p</tt>
  /// \param p A term
  static inline
  void make_minus(term_type& result, const term_type& p)
  {
    pres_system::make_minus(result, p);
  }

  /// \brief Make a conjunction
  /// \param p A term
  /// \param q A term
  /// \return The value <tt>p && q</tt>
  static inline
  term_type and_(const term_type& p, const term_type& q)
  {
    return pres_system::and_(p,q);
  }

  /// \brief Make a conjunction
  /// \param result The value <tt>p && q</tt>
  /// \param p A term
  /// \param q A term
  static inline
  void make_and_(term_type& result, const term_type& p, const term_type& q)
  {
    pres_system::make_and_(result, p, q);
  }

  /// \brief Make a disjunction
  /// \param p A term
  /// \param q A term
  /// \return The value <tt>p || q</tt>
  static inline
  term_type or_(const term_type& p, const term_type& q)
  {
    return pres_system::or_(p,q);
  }

  /// \brief Make a disjunction
  /// \param result The value <tt>p || q</tt>
  /// \param p A term
  /// \param q A term
  static inline
  void make_or_(term_type& result, const term_type& p, const term_type& q)
  {
    pres_system::make_or_(result, p,q);
  }

  template <typename FwdIt>
  static inline
  term_type join_or(FwdIt first, FwdIt last)
  {
    return utilities::detail::join(first, last, or_, false_());
  }

  template <typename FwdIt>
  static inline
  term_type join_and(FwdIt first, FwdIt last)
  {
    return utilities::detail::join(first, last, and_, true_());
  }

  /// \brief Make an implication
  /// \param p A term
  /// \param q A term
  /// \return The value <tt>p => q</tt>
  static inline
  term_type imp(const term_type& p, const term_type& q)
  {
    return pres_system::imp(p,q);
  }

  /// \brief Make an implication
  /// \param result The value <tt>p => q</tt>
  /// \param p A term
  /// \param q A term
  static inline
  void make_imp(term_type& result, const term_type& p, const term_type& q)
  {
    pres_system::make_imp(result, p, q);
  }

  /// \brief Make a plus
  /// \param p A term
  /// \param q A term
  /// \return The value <tt>p + q</tt>
  static inline
  term_type plus(const term_type& p, const term_type& q)
  {
    return pres_system::plus(p,q);
  }

  /// \brief Make a plus
  /// \param result The value <tt>p + q</tt>
  /// \param p A term
  /// \param q A term
  static inline
  void make_plus(term_type& result, const term_type& p, const term_type& q)
  {
    pres_system::make_plus(result, p, q);
  }

  /// \brief Make a const_multiply
  /// \param p A term
  /// \param q A term
  /// \return The value <tt>p * q</tt>
  static inline

  term_type const_multiply(const data::data_expression& p, const term_type& q)
  {
    return pres_system::const_multiply(p, q);
  }

  /// \brief Make a const multiply
  /// \param result The value <tt>p * q</tt>
  /// \param p A term
  /// \param q A term
  static inline
  void make_const_multiply(term_type& result, const data::data_expression& p, const term_type& q)
  {
    pres_system::make_const_multiply(result, p, q);
  }

  /// \brief Make a const multiply alt
  /// \param p A term
  /// \param q A term
  /// \return The value <tt>p * q</tt>
  static inline
  term_type const_multiply_alt(const term_type& p, const data::data_expression& q)
  {
    return pres_system::const_multiply_alt(p,q);
  }

  /// \brief Make a const multiply alt
  /// \param result The value <tt>p * q</tt>
  /// \param p A term
  /// \param q A term
  static inline
  void make_const_multiply_alt(term_type& result, const term_type& p, const data::data_expression& q)
  {
    pres_system::make_const_multiply_alt(result, p, q);
  }

  /// \brief Make a generalized minimum
  /// \param l A sequence of variables
  /// \param p A term
  /// \return The value <tt>inf l.p</tt>
  static inline
  term_type infimum(const variable_sequence_type& l, const term_type& p)
  {
    if (l.empty())
    {
      return p;
    }
    return pres_system::infimum(l, p);
  }

  /// \brief Make a generalized maximum
  /// \param result The value <tt>infimum l.p</tt>
  /// \param l A sequence of variables
  /// \param p A term
  static inline
  void make_infimum(term_type& result, const variable_sequence_type& l, const term_type& p)
  {
    if (l.empty())
    {
      result = p;
      return;
    }
    pres_system::make_infimum(result, l, p);
  }

  /// \brief Make a generalized maximum
  /// \param l A sequence of variables
  /// \param p A term
  /// \return The value <tt>sup l.p</tt>
  static inline
  term_type supremum(const variable_sequence_type& l, const term_type& p)
  {
    if (l.empty())
    {
      return p;
    }
    return pres_system::supremum(l, p);
  }

  /// \brief Make a generalized maximum
  /// \param result The value <tt>supremum l.p</tt>
  /// \param l A sequence of variables
  /// \param p A term
  static inline
  void make_supremum(term_type& result, const variable_sequence_type& l, const term_type& p)
  {
    if (l.empty())
    {
      result = p;
      return;
    }
    pres_system::make_supremum(result, l, p);
  }

  /// \brief Make a generalized sum operator
  /// \param l A sequence of variables
  /// \param p A term
  /// \return The value <tt>sum l.p</tt>
  static inline
  term_type sum(const variable_sequence_type& l, const term_type& p)
  {
    if (l.empty())
    {
      return p;
    }
    return pres_system::sum(l, p);
  }

  /// \brief Make a generalized sum operator
  /// \param result The value <tt>sum l.p</tt>
  /// \param l A sequence of variables
  /// \param p A term
  static inline
  void make_sum(term_type& result, const variable_sequence_type& l, const term_type& p)
  {
    if (l.empty())
    {
      result = p;
      return;
    }
    pres_system::make_sum(result, l, p);
  }

  /// \brief Test for the value true
  /// \param t A term
  /// \return True if it is the value \p true
  static inline
  bool is_true(const term_type& t)
  {
    return data::sort_bool::is_true_function_symbol(t);
  }

  /// \brief Test for the value false
  /// \param t A term
  /// \return True if it is the value \p false
  static inline
  bool is_false(const term_type& t)
  {
    return data::sort_bool::is_false_function_symbol(t);
  }

  /// \brief Test for a minus
  /// \param t A term
  /// \return True if the argument is a minus
  static inline
  bool is_minus(const term_type& t)
  {
    return pres_system::is_minus(t);
  }

  /// \brief Test for a conjunction
  /// \param t A term
  /// \return True if it is a conjunction
  static inline
  bool is_and(const term_type& t)
  {
    return pres_system::is_and(t);
  }

  /// \brief Test for a disjunction
  /// \param t A term
  /// \return True if it is a disjunction
  static inline
  bool is_or(const term_type& t)
  {
    return pres_system::is_or(t);
  }

  /// \brief Test for an implication
  /// \param t A term
  /// \return True if it is an implication
  static inline
  bool is_imp(const term_type& t)
  {
    return pres_system::is_imp(t);
  }

  /// \brief Test for an infimum.
  /// \param t A term
  /// \return True if the argument is an infimum application
  static inline
  bool is_infimum(const term_type& t)
  {
    return pres_system::is_infimum(t);
  }

  /// \brief Test for a max quantification
  /// \param t A term
  /// \return True if t is an maximum  quantification
  static inline
  bool is_supremum(const term_type& t)
  {
    return pres_system::is_supremum(t);
  }

  /// \brief Test for data term
  /// \param t A term
  /// \return True if the term is a data term
  static inline
  bool is_data(const term_type& t)
  {
    return data::is_data_expression(t);
  }

  /// \brief Test for propositional variable instantiation
  /// \param t A term
  /// \return True if the term is a propositional variable instantiation
  static inline
  bool is_prop_var(const term_type& t)
  {
    return pres_system::is_propositional_variable_instantiation(t);
  }

  /// \brief Returns the left argument of a term of type and, or or imp
  /// \param t A term
  /// \return The left argument of the term. Also works for data terms
  static inline
  term_type left(const term_type& t)
  {
    return pres_system::accessors::left(t);
  }

  /// \brief Returns the right argument of a term of type and, or or imp
  /// \param t A term
  /// \return The right argument of the term. Also works for data terms
  static inline
  term_type right(const term_type& t)
  {
    return pres_system::accessors::right(t);
  }

  /// \brief Returns the argument of a term of type not
  /// \param t A term
  static inline
  const term_type& minus_arg(const term_type& t)
  {
    assert(is_pres_minus(t));
    return atermpp::down_cast<term_type>(t[0]);
  }

  /// \brief Returns the quantifier variables of a quantifier expression
  /// \param t A term
  /// \return The requested argument. Doesn't work for data terms
  static inline
  const variable_sequence_type& var(const term_type& t)
  {
    // Forall and exists are not fully supported by the data library
    assert(!data::is_data_expression(t) || (!data::is_abstraction(t)
                                        || (!is_infimum(data::abstraction(t)) && !is_supremum(data::abstraction(t)))));
    assert(is_supremum(t) || is_infimum(t));

    return atermpp::down_cast<variable_sequence_type>(t[0]);
  }

  /// \brief Returns the name of a propositional variable instantiation
  /// \param t A term
  /// \return The name of the propositional variable instantiation
  static inline
  const string_type &name(const term_type& t)
  {
    assert(is_prop_var(t));
    return atermpp::down_cast<string_type>(t[0]);
  }

  /// \brief Returns the parameter list of a propositional variable instantiation
  /// \param t A term
  /// \return The parameter list of the propositional variable instantiation
  static inline
  const data_term_sequence_type &param(const term_type& t)
  {
    assert(is_prop_var(t));
    return atermpp::down_cast<data_term_sequence_type>(t[1]);
  }

  /// \brief Conversion from variable to term
  /// \param v A variable
  /// \return The converted variable
  static inline
  const term_type& variable2term(const variable_type& v)
  {
    return atermpp::down_cast<term_type>(v);
  }

  /// \brief Test if a term is a variable
  /// \param t A term
  /// \return True if the term is a variable
  static inline
  bool is_variable(const term_type& t)
  {
    return data::is_variable(t);
  }

  /// \brief Pretty print function
  /// \param t A term
  /// \return Returns a pretty print representation of the term
  static inline
  std::string pp(const term_type& t)
  {
    return pres_system::pp(t);
  }
};

} // namespace mcrl2::core



namespace std
{

  template <>
  struct hash<mcrl2::pres_system::pres_expression>
  {
    std::size_t operator()(const mcrl2::pres_system::pres_expression& x) const
    {
      return hash<atermpp::aterm>()(x);
    }
  };

  template <>
  struct hash<mcrl2::pres_system::propositional_variable_instantiation>
  {
    std::size_t operator()(const mcrl2::pres_system::propositional_variable_instantiation& x) const
    {
      return hash<atermpp::aterm>()(x);
    }
  };

} // namespace std

#endif // MCRL2_PRES_PRES_EXPRESSION_H
