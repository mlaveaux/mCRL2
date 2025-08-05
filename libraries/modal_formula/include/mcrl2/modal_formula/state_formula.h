// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/modal_formula/state_formula.h
/// \brief Add your file description here.

#ifndef MCRL2_MODAL_FORMULA_STATE_FORMULA_H
#define MCRL2_MODAL_FORMULA_STATE_FORMULA_H

#include "mcrl2/modal_formula/regular_formula.h"

namespace mcrl2::state_formulas
{

//--- start generated classes ---//
/// \\brief A state formula
class state_formula: public atermpp::aterm
{
  public:
    /// \\brief Default constructor X3.
    state_formula()
      : atermpp::aterm(core::detail::default_values::StateFrm)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit state_formula(const atermpp::aterm& term)
      : atermpp::aterm(term)
    {
      assert(core::detail::check_rule_StateFrm(*this));
    }

    /// \\brief Constructor Z6.
    state_formula(const data::data_expression& x)
      : atermpp::aterm(x)
    {}

    /// \\brief Constructor Z6.
    state_formula(const data::untyped_data_parameter& x)
      : atermpp::aterm(x)
    {}

    /// Move semantics
    state_formula(const state_formula&) noexcept = default;
    state_formula(state_formula&&) noexcept = default;
    state_formula& operator=(const state_formula&) noexcept = default;
    state_formula& operator=(state_formula&&) noexcept = default;
//--- start user section state_formula ---//
    /// \brief Returns true if the formula is timed.
    /// \return True if the formula is timed.
    bool has_time() const;
//--- end user section state_formula ---//
};

/// \\brief list of state_formulas
using state_formula_list = atermpp::term_list<state_formula>;

/// \\brief vector of state_formulas
using state_formula_vector = std::vector<state_formula>;

// prototypes
inline bool is_true(const atermpp::aterm& x);
inline bool is_false(const atermpp::aterm& x);
inline bool is_not(const atermpp::aterm& x);
inline bool is_minus(const atermpp::aterm& x);
inline bool is_and(const atermpp::aterm& x);
inline bool is_or(const atermpp::aterm& x);
inline bool is_imp(const atermpp::aterm& x);
inline bool is_plus(const atermpp::aterm& x);
inline bool is_const_multiply(const atermpp::aterm& x);
inline bool is_const_multiply_alt(const atermpp::aterm& x);
inline bool is_forall(const atermpp::aterm& x);
inline bool is_exists(const atermpp::aterm& x);
inline bool is_infimum(const atermpp::aterm& x);
inline bool is_supremum(const atermpp::aterm& x);
inline bool is_sum(const atermpp::aterm& x);
inline bool is_must(const atermpp::aterm& x);
inline bool is_may(const atermpp::aterm& x);
inline bool is_yaled(const atermpp::aterm& x);
inline bool is_yaled_timed(const atermpp::aterm& x);
inline bool is_delay(const atermpp::aterm& x);
inline bool is_delay_timed(const atermpp::aterm& x);
inline bool is_variable(const atermpp::aterm& x);
inline bool is_nu(const atermpp::aterm& x);
inline bool is_mu(const atermpp::aterm& x);

/// \\brief Test for a state_formula expression
/// \\param x A term
/// \\return True if \\a x is a state_formula expression
inline
bool is_state_formula(const atermpp::aterm& x)
{
  return data::is_data_expression(x) ||
         data::is_untyped_data_parameter(x) ||
         state_formulas::is_true(x) ||
         state_formulas::is_false(x) ||
         state_formulas::is_not(x) ||
         state_formulas::is_minus(x) ||
         state_formulas::is_and(x) ||
         state_formulas::is_or(x) ||
         state_formulas::is_imp(x) ||
         state_formulas::is_plus(x) ||
         state_formulas::is_const_multiply(x) ||
         state_formulas::is_const_multiply_alt(x) ||
         state_formulas::is_forall(x) ||
         state_formulas::is_exists(x) ||
         state_formulas::is_infimum(x) ||
         state_formulas::is_supremum(x) ||
         state_formulas::is_sum(x) ||
         state_formulas::is_must(x) ||
         state_formulas::is_may(x) ||
         state_formulas::is_yaled(x) ||
         state_formulas::is_yaled_timed(x) ||
         state_formulas::is_delay(x) ||
         state_formulas::is_delay_timed(x) ||
         state_formulas::is_variable(x) ||
         state_formulas::is_nu(x) ||
         state_formulas::is_mu(x);
}

// prototype declaration
std::string pp(const state_formula& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const state_formula& x)
{
  return out << state_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(state_formula& t1, state_formula& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The value true for state formulas
class true_: public state_formula
{
  public:
    /// \\brief Default constructor X3.
    true_()
      : state_formula(core::detail::default_values::StateTrue)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit true_(const atermpp::aterm& term)
      : state_formula(term)
    {
      assert(core::detail::check_term_StateTrue(*this));
    }

    /// Move semantics
    true_(const true_&) noexcept = default;
    true_(true_&&) noexcept = default;
    true_& operator=(const true_&) noexcept = default;
    true_& operator=(true_&&) noexcept = default;
};

/// \\brief Test for a true expression
/// \\param x A term
/// \\return True if \\a x is a true expression
inline
bool is_true(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::StateTrue;
}

// prototype declaration
std::string pp(const true_& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const true_& x)
{
  return out << state_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(true_& t1, true_& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The value false for state formulas
class false_: public state_formula
{
  public:
    /// \\brief Default constructor X3.
    false_()
      : state_formula(core::detail::default_values::StateFalse)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit false_(const atermpp::aterm& term)
      : state_formula(term)
    {
      assert(core::detail::check_term_StateFalse(*this));
    }

    /// Move semantics
    false_(const false_&) noexcept = default;
    false_(false_&&) noexcept = default;
    false_& operator=(const false_&) noexcept = default;
    false_& operator=(false_&&) noexcept = default;
};

/// \\brief Test for a false expression
/// \\param x A term
/// \\return True if \\a x is a false expression
inline
bool is_false(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::StateFalse;
}

// prototype declaration
std::string pp(const false_& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const false_& x)
{
  return out << state_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(false_& t1, false_& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The not operator for state formulas
class not_: public state_formula
{
  public:
    /// \\brief Default constructor X3.
    not_()
      : state_formula(core::detail::default_values::StateNot)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit not_(const atermpp::aterm& term)
      : state_formula(term)
    {
      assert(core::detail::check_term_StateNot(*this));
    }

    /// \\brief Constructor Z14.
    explicit not_(const state_formula& operand)
      : state_formula(atermpp::aterm(core::detail::function_symbol_StateNot(), operand))
    {}

    /// Move semantics
    not_(const not_&) noexcept = default;
    not_(not_&&) noexcept = default;
    not_& operator=(const not_&) noexcept = default;
    not_& operator=(not_&&) noexcept = default;

    const state_formula& operand() const
    {
      return atermpp::down_cast<state_formula>((*this)[0]);
    }
};

/// \\brief Make_not_ constructs a new term into a given address.
/// \\ \param t The reference into which the new not_ is constructed. 
template <class... ARGUMENTS>
inline void make_not_(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_StateNot(), args...);
}

/// \\brief Test for a not expression
/// \\param x A term
/// \\return True if \\a x is a not expression
inline
bool is_not(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::StateNot;
}

// prototype declaration
std::string pp(const not_& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const not_& x)
{
  return out << state_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(not_& t1, not_& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The minus operator for state formulas
class minus: public state_formula
{
  public:
    /// \\brief Default constructor X3.
    minus()
      : state_formula(core::detail::default_values::StateMinus)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit minus(const atermpp::aterm& term)
      : state_formula(term)
    {
      assert(core::detail::check_term_StateMinus(*this));
    }

    /// \\brief Constructor Z14.
    explicit minus(const state_formula& operand)
      : state_formula(atermpp::aterm(core::detail::function_symbol_StateMinus(), operand))
    {}

    /// Move semantics
    minus(const minus&) noexcept = default;
    minus(minus&&) noexcept = default;
    minus& operator=(const minus&) noexcept = default;
    minus& operator=(minus&&) noexcept = default;

    const state_formula& operand() const
    {
      return atermpp::down_cast<state_formula>((*this)[0]);
    }
};

/// \\brief Make_minus constructs a new term into a given address.
/// \\ \param t The reference into which the new minus is constructed. 
template <class... ARGUMENTS>
inline void make_minus(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_StateMinus(), args...);
}

/// \\brief Test for a minus expression
/// \\param x A term
/// \\return True if \\a x is a minus expression
inline
bool is_minus(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::StateMinus;
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
  return out << state_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(minus& t1, minus& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The and operator for state formulas
class and_: public state_formula
{
  public:
    /// \\brief Default constructor X3.
    and_()
      : state_formula(core::detail::default_values::StateAnd)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit and_(const atermpp::aterm& term)
      : state_formula(term)
    {
      assert(core::detail::check_term_StateAnd(*this));
    }

    /// \\brief Constructor Z14.
    and_(const state_formula& left, const state_formula& right)
      : state_formula(atermpp::aterm(core::detail::function_symbol_StateAnd(), left, right))
    {}

    /// Move semantics
    and_(const and_&) noexcept = default;
    and_(and_&&) noexcept = default;
    and_& operator=(const and_&) noexcept = default;
    and_& operator=(and_&&) noexcept = default;

    const state_formula& left() const
    {
      return atermpp::down_cast<state_formula>((*this)[0]);
    }

    const state_formula& right() const
    {
      return atermpp::down_cast<state_formula>((*this)[1]);
    }
};

/// \\brief Make_and_ constructs a new term into a given address.
/// \\ \param t The reference into which the new and_ is constructed. 
template <class... ARGUMENTS>
inline void make_and_(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_StateAnd(), args...);
}

/// \\brief Test for a and expression
/// \\param x A term
/// \\return True if \\a x is a and expression
inline
bool is_and(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::StateAnd;
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
  return out << state_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(and_& t1, and_& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The or operator for state formulas
class or_: public state_formula
{
  public:
    /// \\brief Default constructor X3.
    or_()
      : state_formula(core::detail::default_values::StateOr)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit or_(const atermpp::aterm& term)
      : state_formula(term)
    {
      assert(core::detail::check_term_StateOr(*this));
    }

    /// \\brief Constructor Z14.
    or_(const state_formula& left, const state_formula& right)
      : state_formula(atermpp::aterm(core::detail::function_symbol_StateOr(), left, right))
    {}

    /// Move semantics
    or_(const or_&) noexcept = default;
    or_(or_&&) noexcept = default;
    or_& operator=(const or_&) noexcept = default;
    or_& operator=(or_&&) noexcept = default;

    const state_formula& left() const
    {
      return atermpp::down_cast<state_formula>((*this)[0]);
    }

    const state_formula& right() const
    {
      return atermpp::down_cast<state_formula>((*this)[1]);
    }
};

/// \\brief Make_or_ constructs a new term into a given address.
/// \\ \param t The reference into which the new or_ is constructed. 
template <class... ARGUMENTS>
inline void make_or_(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_StateOr(), args...);
}

/// \\brief Test for a or expression
/// \\param x A term
/// \\return True if \\a x is a or expression
inline
bool is_or(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::StateOr;
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
  return out << state_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(or_& t1, or_& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The implication operator for state formulas
class imp: public state_formula
{
  public:
    /// \\brief Default constructor X3.
    imp()
      : state_formula(core::detail::default_values::StateImp)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit imp(const atermpp::aterm& term)
      : state_formula(term)
    {
      assert(core::detail::check_term_StateImp(*this));
    }

    /// \\brief Constructor Z14.
    imp(const state_formula& left, const state_formula& right)
      : state_formula(atermpp::aterm(core::detail::function_symbol_StateImp(), left, right))
    {}

    /// Move semantics
    imp(const imp&) noexcept = default;
    imp(imp&&) noexcept = default;
    imp& operator=(const imp&) noexcept = default;
    imp& operator=(imp&&) noexcept = default;

    const state_formula& left() const
    {
      return atermpp::down_cast<state_formula>((*this)[0]);
    }

    const state_formula& right() const
    {
      return atermpp::down_cast<state_formula>((*this)[1]);
    }
};

/// \\brief Make_imp constructs a new term into a given address.
/// \\ \param t The reference into which the new imp is constructed. 
template <class... ARGUMENTS>
inline void make_imp(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_StateImp(), args...);
}

/// \\brief Test for a imp expression
/// \\param x A term
/// \\return True if \\a x is a imp expression
inline
bool is_imp(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::StateImp;
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
  return out << state_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(imp& t1, imp& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The plus operator for state formulas with values
class plus: public state_formula
{
  public:
    /// \\brief Default constructor X3.
    plus()
      : state_formula(core::detail::default_values::StatePlus)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit plus(const atermpp::aterm& term)
      : state_formula(term)
    {
      assert(core::detail::check_term_StatePlus(*this));
    }

    /// \\brief Constructor Z14.
    plus(const state_formula& left, const state_formula& right)
      : state_formula(atermpp::aterm(core::detail::function_symbol_StatePlus(), left, right))
    {}

    /// Move semantics
    plus(const plus&) noexcept = default;
    plus(plus&&) noexcept = default;
    plus& operator=(const plus&) noexcept = default;
    plus& operator=(plus&&) noexcept = default;

    const state_formula& left() const
    {
      return atermpp::down_cast<state_formula>((*this)[0]);
    }

    const state_formula& right() const
    {
      return atermpp::down_cast<state_formula>((*this)[1]);
    }
};

/// \\brief Make_plus constructs a new term into a given address.
/// \\ \param t The reference into which the new plus is constructed. 
template <class... ARGUMENTS>
inline void make_plus(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_StatePlus(), args...);
}

/// \\brief Test for a plus expression
/// \\param x A term
/// \\return True if \\a x is a plus expression
inline
bool is_plus(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::StatePlus;
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
  return out << state_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(plus& t1, plus& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The multiply operator for state formulas with values
class const_multiply: public state_formula
{
  public:
    /// \\brief Default constructor X3.
    const_multiply()
      : state_formula(core::detail::default_values::StateConstantMultiply)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit const_multiply(const atermpp::aterm& term)
      : state_formula(term)
    {
      assert(core::detail::check_term_StateConstantMultiply(*this));
    }

    /// \\brief Constructor Z14.
    const_multiply(const data::data_expression& left, const state_formula& right)
      : state_formula(atermpp::aterm(core::detail::function_symbol_StateConstantMultiply(), left, right))
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

    const state_formula& right() const
    {
      return atermpp::down_cast<state_formula>((*this)[1]);
    }
};

/// \\brief Make_const_multiply constructs a new term into a given address.
/// \\ \param t The reference into which the new const_multiply is constructed. 
template <class... ARGUMENTS>
inline void make_const_multiply(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_StateConstantMultiply(), args...);
}

/// \\brief Test for a const_multiply expression
/// \\param x A term
/// \\return True if \\a x is a const_multiply expression
inline
bool is_const_multiply(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::StateConstantMultiply;
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
  return out << state_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(const_multiply& t1, const_multiply& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The multiply operator for state formulas with values
class const_multiply_alt: public state_formula
{
  public:
    /// \\brief Default constructor X3.
    const_multiply_alt()
      : state_formula(core::detail::default_values::StateConstantMultiplyAlt)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit const_multiply_alt(const atermpp::aterm& term)
      : state_formula(term)
    {
      assert(core::detail::check_term_StateConstantMultiplyAlt(*this));
    }

    /// \\brief Constructor Z14.
    const_multiply_alt(const state_formula& left, const data::data_expression& right)
      : state_formula(atermpp::aterm(core::detail::function_symbol_StateConstantMultiplyAlt(), left, right))
    {}

    /// Move semantics
    const_multiply_alt(const const_multiply_alt&) noexcept = default;
    const_multiply_alt(const_multiply_alt&&) noexcept = default;
    const_multiply_alt& operator=(const const_multiply_alt&) noexcept = default;
    const_multiply_alt& operator=(const_multiply_alt&&) noexcept = default;

    const state_formula& left() const
    {
      return atermpp::down_cast<state_formula>((*this)[0]);
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
  atermpp::make_term_appl(t, core::detail::function_symbol_StateConstantMultiplyAlt(), args...);
}

/// \\brief Test for a const_multiply_alt expression
/// \\param x A term
/// \\return True if \\a x is a const_multiply_alt expression
inline
bool is_const_multiply_alt(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::StateConstantMultiplyAlt;
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
  return out << state_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(const_multiply_alt& t1, const_multiply_alt& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The universal quantification operator for state formulas
class forall: public state_formula
{
  public:
    /// \\brief Default constructor X3.
    forall()
      : state_formula(core::detail::default_values::StateForall)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit forall(const atermpp::aterm& term)
      : state_formula(term)
    {
      assert(core::detail::check_term_StateForall(*this));
    }

    /// \\brief Constructor Z14.
    forall(const data::variable_list& variables, const state_formula& body)
      : state_formula(atermpp::aterm(core::detail::function_symbol_StateForall(), variables, body))
    {}

    /// Move semantics
    forall(const forall&) noexcept = default;
    forall(forall&&) noexcept = default;
    forall& operator=(const forall&) noexcept = default;
    forall& operator=(forall&&) noexcept = default;

    const data::variable_list& variables() const
    {
      return atermpp::down_cast<data::variable_list>((*this)[0]);
    }

    const state_formula& body() const
    {
      return atermpp::down_cast<state_formula>((*this)[1]);
    }
};

/// \\brief Make_forall constructs a new term into a given address.
/// \\ \param t The reference into which the new forall is constructed. 
template <class... ARGUMENTS>
inline void make_forall(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_StateForall(), args...);
}

/// \\brief Test for a forall expression
/// \\param x A term
/// \\return True if \\a x is a forall expression
inline
bool is_forall(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::StateForall;
}

// prototype declaration
std::string pp(const forall& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const forall& x)
{
  return out << state_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(forall& t1, forall& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The existential quantification operator for state formulas
class exists: public state_formula
{
  public:
    /// \\brief Default constructor X3.
    exists()
      : state_formula(core::detail::default_values::StateExists)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit exists(const atermpp::aterm& term)
      : state_formula(term)
    {
      assert(core::detail::check_term_StateExists(*this));
    }

    /// \\brief Constructor Z14.
    exists(const data::variable_list& variables, const state_formula& body)
      : state_formula(atermpp::aterm(core::detail::function_symbol_StateExists(), variables, body))
    {}

    /// Move semantics
    exists(const exists&) noexcept = default;
    exists(exists&&) noexcept = default;
    exists& operator=(const exists&) noexcept = default;
    exists& operator=(exists&&) noexcept = default;

    const data::variable_list& variables() const
    {
      return atermpp::down_cast<data::variable_list>((*this)[0]);
    }

    const state_formula& body() const
    {
      return atermpp::down_cast<state_formula>((*this)[1]);
    }
};

/// \\brief Make_exists constructs a new term into a given address.
/// \\ \param t The reference into which the new exists is constructed. 
template <class... ARGUMENTS>
inline void make_exists(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_StateExists(), args...);
}

/// \\brief Test for a exists expression
/// \\param x A term
/// \\return True if \\a x is a exists expression
inline
bool is_exists(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::StateExists;
}

// prototype declaration
std::string pp(const exists& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const exists& x)
{
  return out << state_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(exists& t1, exists& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The infimum over a data type for state formulas
class infimum: public state_formula
{
  public:
    /// \\brief Default constructor X3.
    infimum()
      : state_formula(core::detail::default_values::StateInfimum)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit infimum(const atermpp::aterm& term)
      : state_formula(term)
    {
      assert(core::detail::check_term_StateInfimum(*this));
    }

    /// \\brief Constructor Z14.
    infimum(const data::variable_list& variables, const state_formula& body)
      : state_formula(atermpp::aterm(core::detail::function_symbol_StateInfimum(), variables, body))
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

    const state_formula& body() const
    {
      return atermpp::down_cast<state_formula>((*this)[1]);
    }
};

/// \\brief Make_infimum constructs a new term into a given address.
/// \\ \param t The reference into which the new infimum is constructed. 
template <class... ARGUMENTS>
inline void make_infimum(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_StateInfimum(), args...);
}

/// \\brief Test for a infimum expression
/// \\param x A term
/// \\return True if \\a x is a infimum expression
inline
bool is_infimum(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::StateInfimum;
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
  return out << state_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(infimum& t1, infimum& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The supremum over a data type for state formulas
class supremum: public state_formula
{
  public:
    /// \\brief Default constructor X3.
    supremum()
      : state_formula(core::detail::default_values::StateSupremum)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit supremum(const atermpp::aterm& term)
      : state_formula(term)
    {
      assert(core::detail::check_term_StateSupremum(*this));
    }

    /// \\brief Constructor Z14.
    supremum(const data::variable_list& variables, const state_formula& body)
      : state_formula(atermpp::aterm(core::detail::function_symbol_StateSupremum(), variables, body))
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

    const state_formula& body() const
    {
      return atermpp::down_cast<state_formula>((*this)[1]);
    }
};

/// \\brief Make_supremum constructs a new term into a given address.
/// \\ \param t The reference into which the new supremum is constructed. 
template <class... ARGUMENTS>
inline void make_supremum(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_StateSupremum(), args...);
}

/// \\brief Test for a supremum expression
/// \\param x A term
/// \\return True if \\a x is a supremum expression
inline
bool is_supremum(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::StateSupremum;
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
  return out << state_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(supremum& t1, supremum& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The sum over a data type for state formulas
class sum: public state_formula
{
  public:
    /// \\brief Default constructor X3.
    sum()
      : state_formula(core::detail::default_values::StateSum)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit sum(const atermpp::aterm& term)
      : state_formula(term)
    {
      assert(core::detail::check_term_StateSum(*this));
    }

    /// \\brief Constructor Z14.
    sum(const data::variable_list& variables, const state_formula& body)
      : state_formula(atermpp::aterm(core::detail::function_symbol_StateSum(), variables, body))
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

    const state_formula& body() const
    {
      return atermpp::down_cast<state_formula>((*this)[1]);
    }
};

/// \\brief Make_sum constructs a new term into a given address.
/// \\ \param t The reference into which the new sum is constructed. 
template <class... ARGUMENTS>
inline void make_sum(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_StateSum(), args...);
}

/// \\brief Test for a sum expression
/// \\param x A term
/// \\return True if \\a x is a sum expression
inline
bool is_sum(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::StateSum;
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
  return out << state_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(sum& t1, sum& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The must operator for state formulas
class must: public state_formula
{
  public:
    /// \\brief Default constructor X3.
    must()
      : state_formula(core::detail::default_values::StateMust)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit must(const atermpp::aterm& term)
      : state_formula(term)
    {
      assert(core::detail::check_term_StateMust(*this));
    }

    /// \\brief Constructor Z14.
    must(const regular_formulas::regular_formula& formula, const state_formula& operand)
      : state_formula(atermpp::aterm(core::detail::function_symbol_StateMust(), formula, operand))
    {}

    /// Move semantics
    must(const must&) noexcept = default;
    must(must&&) noexcept = default;
    must& operator=(const must&) noexcept = default;
    must& operator=(must&&) noexcept = default;

    const regular_formulas::regular_formula& formula() const
    {
      return atermpp::down_cast<regular_formulas::regular_formula>((*this)[0]);
    }

    const state_formula& operand() const
    {
      return atermpp::down_cast<state_formula>((*this)[1]);
    }
};

/// \\brief Make_must constructs a new term into a given address.
/// \\ \param t The reference into which the new must is constructed. 
template <class... ARGUMENTS>
inline void make_must(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_StateMust(), args...);
}

/// \\brief Test for a must expression
/// \\param x A term
/// \\return True if \\a x is a must expression
inline
bool is_must(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::StateMust;
}

// prototype declaration
std::string pp(const must& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const must& x)
{
  return out << state_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(must& t1, must& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The may operator for state formulas
class may: public state_formula
{
  public:
    /// \\brief Default constructor X3.
    may()
      : state_formula(core::detail::default_values::StateMay)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit may(const atermpp::aterm& term)
      : state_formula(term)
    {
      assert(core::detail::check_term_StateMay(*this));
    }

    /// \\brief Constructor Z14.
    may(const regular_formulas::regular_formula& formula, const state_formula& operand)
      : state_formula(atermpp::aterm(core::detail::function_symbol_StateMay(), formula, operand))
    {}

    /// Move semantics
    may(const may&) noexcept = default;
    may(may&&) noexcept = default;
    may& operator=(const may&) noexcept = default;
    may& operator=(may&&) noexcept = default;

    const regular_formulas::regular_formula& formula() const
    {
      return atermpp::down_cast<regular_formulas::regular_formula>((*this)[0]);
    }

    const state_formula& operand() const
    {
      return atermpp::down_cast<state_formula>((*this)[1]);
    }
};

/// \\brief Make_may constructs a new term into a given address.
/// \\ \param t The reference into which the new may is constructed. 
template <class... ARGUMENTS>
inline void make_may(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_StateMay(), args...);
}

/// \\brief Test for a may expression
/// \\param x A term
/// \\return True if \\a x is a may expression
inline
bool is_may(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::StateMay;
}

// prototype declaration
std::string pp(const may& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const may& x)
{
  return out << state_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(may& t1, may& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The yaled operator for state formulas
class yaled: public state_formula
{
  public:
    /// \\brief Default constructor X3.
    yaled()
      : state_formula(core::detail::default_values::StateYaled)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit yaled(const atermpp::aterm& term)
      : state_formula(term)
    {
      assert(core::detail::check_term_StateYaled(*this));
    }

    /// Move semantics
    yaled(const yaled&) noexcept = default;
    yaled(yaled&&) noexcept = default;
    yaled& operator=(const yaled&) noexcept = default;
    yaled& operator=(yaled&&) noexcept = default;
};

/// \\brief Test for a yaled expression
/// \\param x A term
/// \\return True if \\a x is a yaled expression
inline
bool is_yaled(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::StateYaled;
}

// prototype declaration
std::string pp(const yaled& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const yaled& x)
{
  return out << state_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(yaled& t1, yaled& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The timed yaled operator for state formulas
class yaled_timed: public state_formula
{
  public:
    /// \\brief Default constructor X3.
    yaled_timed()
      : state_formula(core::detail::default_values::StateYaledTimed)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit yaled_timed(const atermpp::aterm& term)
      : state_formula(term)
    {
      assert(core::detail::check_term_StateYaledTimed(*this));
    }

    /// \\brief Constructor Z14.
    explicit yaled_timed(const data::data_expression& time_stamp)
      : state_formula(atermpp::aterm(core::detail::function_symbol_StateYaledTimed(), time_stamp))
    {}

    /// Move semantics
    yaled_timed(const yaled_timed&) noexcept = default;
    yaled_timed(yaled_timed&&) noexcept = default;
    yaled_timed& operator=(const yaled_timed&) noexcept = default;
    yaled_timed& operator=(yaled_timed&&) noexcept = default;

    const data::data_expression& time_stamp() const
    {
      return atermpp::down_cast<data::data_expression>((*this)[0]);
    }
};

/// \\brief Make_yaled_timed constructs a new term into a given address.
/// \\ \param t The reference into which the new yaled_timed is constructed. 
template <class... ARGUMENTS>
inline void make_yaled_timed(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_StateYaledTimed(), args...);
}

/// \\brief Test for a yaled_timed expression
/// \\param x A term
/// \\return True if \\a x is a yaled_timed expression
inline
bool is_yaled_timed(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::StateYaledTimed;
}

// prototype declaration
std::string pp(const yaled_timed& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const yaled_timed& x)
{
  return out << state_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(yaled_timed& t1, yaled_timed& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The delay operator for state formulas
class delay: public state_formula
{
  public:
    /// \\brief Default constructor X3.
    delay()
      : state_formula(core::detail::default_values::StateDelay)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit delay(const atermpp::aterm& term)
      : state_formula(term)
    {
      assert(core::detail::check_term_StateDelay(*this));
    }

    /// Move semantics
    delay(const delay&) noexcept = default;
    delay(delay&&) noexcept = default;
    delay& operator=(const delay&) noexcept = default;
    delay& operator=(delay&&) noexcept = default;
};

/// \\brief Test for a delay expression
/// \\param x A term
/// \\return True if \\a x is a delay expression
inline
bool is_delay(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::StateDelay;
}

// prototype declaration
std::string pp(const delay& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const delay& x)
{
  return out << state_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(delay& t1, delay& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The timed delay operator for state formulas
class delay_timed: public state_formula
{
  public:
    /// \\brief Default constructor X3.
    delay_timed()
      : state_formula(core::detail::default_values::StateDelayTimed)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit delay_timed(const atermpp::aterm& term)
      : state_formula(term)
    {
      assert(core::detail::check_term_StateDelayTimed(*this));
    }

    /// \\brief Constructor Z14.
    explicit delay_timed(const data::data_expression& time_stamp)
      : state_formula(atermpp::aterm(core::detail::function_symbol_StateDelayTimed(), time_stamp))
    {}

    /// Move semantics
    delay_timed(const delay_timed&) noexcept = default;
    delay_timed(delay_timed&&) noexcept = default;
    delay_timed& operator=(const delay_timed&) noexcept = default;
    delay_timed& operator=(delay_timed&&) noexcept = default;

    const data::data_expression& time_stamp() const
    {
      return atermpp::down_cast<data::data_expression>((*this)[0]);
    }
};

/// \\brief Make_delay_timed constructs a new term into a given address.
/// \\ \param t The reference into which the new delay_timed is constructed. 
template <class... ARGUMENTS>
inline void make_delay_timed(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_StateDelayTimed(), args...);
}

/// \\brief Test for a delay_timed expression
/// \\param x A term
/// \\return True if \\a x is a delay_timed expression
inline
bool is_delay_timed(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::StateDelayTimed;
}

// prototype declaration
std::string pp(const delay_timed& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const delay_timed& x)
{
  return out << state_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(delay_timed& t1, delay_timed& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The state formula variable
class variable: public state_formula
{
  public:
    /// \\brief Default constructor X3.
    variable()
      : state_formula(core::detail::default_values::StateVar)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit variable(const atermpp::aterm& term)
      : state_formula(term)
    {
      assert(core::detail::check_term_StateVar(*this));
    }

    /// \\brief Constructor Z14.
    variable(const core::identifier_string& name, const data::data_expression_list& arguments)
      : state_formula(atermpp::aterm(core::detail::function_symbol_StateVar(), name, arguments))
    {}

    /// \\brief Constructor Z2.
    variable(const std::string& name, const data::data_expression_list& arguments)
      : state_formula(atermpp::aterm(core::detail::function_symbol_StateVar(), core::identifier_string(name), arguments))
    {}

    /// Move semantics
    variable(const variable&) noexcept = default;
    variable(variable&&) noexcept = default;
    variable& operator=(const variable&) noexcept = default;
    variable& operator=(variable&&) noexcept = default;

    const core::identifier_string& name() const
    {
      return atermpp::down_cast<core::identifier_string>((*this)[0]);
    }

    const data::data_expression_list& arguments() const
    {
      return atermpp::down_cast<data::data_expression_list>((*this)[1]);
    }
};

/// \\brief Make_variable constructs a new term into a given address.
/// \\ \param t The reference into which the new variable is constructed. 
template <class... ARGUMENTS>
inline void make_variable(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_StateVar(), args...);
}

/// \\brief Test for a variable expression
/// \\param x A term
/// \\return True if \\a x is a variable expression
inline
bool is_variable(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::StateVar;
}

// prototype declaration
std::string pp(const variable& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const variable& x)
{
  return out << state_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(variable& t1, variable& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The nu operator for state formulas
class nu: public state_formula
{
  public:
    /// \\brief Default constructor X3.
    nu()
      : state_formula(core::detail::default_values::StateNu)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit nu(const atermpp::aterm& term)
      : state_formula(term)
    {
      assert(core::detail::check_term_StateNu(*this));
    }

    /// \\brief Constructor Z14.
    nu(const core::identifier_string& name, const data::assignment_list& assignments, const state_formula& operand)
      : state_formula(atermpp::aterm(core::detail::function_symbol_StateNu(), name, assignments, operand))
    {}

    /// \\brief Constructor Z2.
    nu(const std::string& name, const data::assignment_list& assignments, const state_formula& operand)
      : state_formula(atermpp::aterm(core::detail::function_symbol_StateNu(), core::identifier_string(name), assignments, operand))
    {}

    /// Move semantics
    nu(const nu&) noexcept = default;
    nu(nu&&) noexcept = default;
    nu& operator=(const nu&) noexcept = default;
    nu& operator=(nu&&) noexcept = default;

    const core::identifier_string& name() const
    {
      return atermpp::down_cast<core::identifier_string>((*this)[0]);
    }

    const data::assignment_list& assignments() const
    {
      return atermpp::down_cast<data::assignment_list>((*this)[1]);
    }

    const state_formula& operand() const
    {
      return atermpp::down_cast<state_formula>((*this)[2]);
    }
};

/// \\brief Make_nu constructs a new term into a given address.
/// \\ \param t The reference into which the new nu is constructed. 
template <class... ARGUMENTS>
inline void make_nu(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_StateNu(), args...);
}

/// \\brief Test for a nu expression
/// \\param x A term
/// \\return True if \\a x is a nu expression
inline
bool is_nu(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::StateNu;
}

// prototype declaration
std::string pp(const nu& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const nu& x)
{
  return out << state_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(nu& t1, nu& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The mu operator for state formulas
class mu: public state_formula
{
  public:
    /// \\brief Default constructor X3.
    mu()
      : state_formula(core::detail::default_values::StateMu)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit mu(const atermpp::aterm& term)
      : state_formula(term)
    {
      assert(core::detail::check_term_StateMu(*this));
    }

    /// \\brief Constructor Z14.
    mu(const core::identifier_string& name, const data::assignment_list& assignments, const state_formula& operand)
      : state_formula(atermpp::aterm(core::detail::function_symbol_StateMu(), name, assignments, operand))
    {}

    /// \\brief Constructor Z2.
    mu(const std::string& name, const data::assignment_list& assignments, const state_formula& operand)
      : state_formula(atermpp::aterm(core::detail::function_symbol_StateMu(), core::identifier_string(name), assignments, operand))
    {}

    /// Move semantics
    mu(const mu&) noexcept = default;
    mu(mu&&) noexcept = default;
    mu& operator=(const mu&) noexcept = default;
    mu& operator=(mu&&) noexcept = default;

    const core::identifier_string& name() const
    {
      return atermpp::down_cast<core::identifier_string>((*this)[0]);
    }

    const data::assignment_list& assignments() const
    {
      return atermpp::down_cast<data::assignment_list>((*this)[1]);
    }

    const state_formula& operand() const
    {
      return atermpp::down_cast<state_formula>((*this)[2]);
    }
};

/// \\brief Make_mu constructs a new term into a given address.
/// \\ \param t The reference into which the new mu is constructed. 
template <class... ARGUMENTS>
inline void make_mu(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_StateMu(), args...);
}

/// \\brief Test for a mu expression
/// \\param x A term
/// \\return True if \\a x is a mu expression
inline
bool is_mu(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::StateMu;
}

// prototype declaration
std::string pp(const mu& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const mu& x)
{
  return out << state_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(mu& t1, mu& t2) noexcept
{
  t1.swap(t2);
}
//--- end generated classes ---//

namespace algorithms {
    bool is_timed(const state_formula& x);
} // namespace algorithms

/// \brief Returns true if the formula is timed.
/// \return True if the formula is timed.
inline
bool state_formula::has_time() const
{
  return algorithms::is_timed(*this);
}

// template function overloads
state_formula normalize_sorts(const state_formula& x, const data::sort_specification& sortspec);
state_formulas::state_formula translate_user_notation(const state_formulas::state_formula& x);
std::set<data::sort_expression> find_sort_expressions(const state_formulas::state_formula& x);
std::set<data::variable> find_all_variables(const state_formulas::state_formula& x);
std::set<data::variable> find_free_variables(const state_formulas::state_formula& x);
std::set<core::identifier_string> find_identifiers(const state_formulas::state_formula& x);
std::set<process::action_label> find_action_labels(const state_formulas::state_formula& x);
bool find_nil(const state_formulas::state_formula& x);

} // namespace mcrl2::state_formulas



#endif // MCRL2_MODAL_FORMULA_STATE_FORMULA_H
