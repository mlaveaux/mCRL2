// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/modal_formula/action_formula.h
/// \brief Add your file description here.

#ifndef MCRL2_MODAL_FORMULA_ACTION_FORMULA_H
#define MCRL2_MODAL_FORMULA_ACTION_FORMULA_H

#include "mcrl2/lps/multi_action.h"
#include "mcrl2/process/untyped_multi_action.h"

namespace mcrl2::action_formulas
{

//--- start generated classes ---//
/// \\brief An action formula
class action_formula: public atermpp::aterm
{
  public:
    /// \\brief Default constructor X3.
    action_formula()
      : atermpp::aterm(core::detail::default_values::ActFrm)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit action_formula(const atermpp::aterm& term)
      : atermpp::aterm(term)
    {
      assert(core::detail::check_rule_ActFrm(*this));
    }

    /// \\brief Constructor Z6.
    action_formula(const data::data_expression& x)
      : atermpp::aterm(x)
    {}

    /// \\brief Constructor Z6.
    action_formula(const data::untyped_data_parameter& x)
      : atermpp::aterm(x)
    {}

    /// \\brief Constructor Z6.
    action_formula(const process::untyped_multi_action& x)
      : atermpp::aterm(x)
    {}

    /// Move semantics
    action_formula(const action_formula&) noexcept = default;
    action_formula(action_formula&&) noexcept = default;
    action_formula& operator=(const action_formula&) noexcept = default;
    action_formula& operator=(action_formula&&) noexcept = default;
};

/// \\brief list of action_formulas
using action_formula_list = atermpp::term_list<action_formula>;

/// \\brief vector of action_formulas
using action_formula_vector = std::vector<action_formula>;

// prototypes
inline bool is_true(const atermpp::aterm& x);
inline bool is_false(const atermpp::aterm& x);
inline bool is_not(const atermpp::aterm& x);
inline bool is_and(const atermpp::aterm& x);
inline bool is_or(const atermpp::aterm& x);
inline bool is_imp(const atermpp::aterm& x);
inline bool is_forall(const atermpp::aterm& x);
inline bool is_exists(const atermpp::aterm& x);
inline bool is_at(const atermpp::aterm& x);
inline bool is_multi_action(const atermpp::aterm& x);

/// \\brief Test for a action_formula expression
/// \\param x A term
/// \\return True if \\a x is a action_formula expression
inline
bool is_action_formula(const atermpp::aterm& x)
{
  return data::is_data_expression(x) ||
         data::is_untyped_data_parameter(x) ||
         action_formulas::is_true(x) ||
         action_formulas::is_false(x) ||
         action_formulas::is_not(x) ||
         action_formulas::is_and(x) ||
         action_formulas::is_or(x) ||
         action_formulas::is_imp(x) ||
         action_formulas::is_forall(x) ||
         action_formulas::is_exists(x) ||
         action_formulas::is_at(x) ||
         action_formulas::is_multi_action(x) ||
         process::is_untyped_multi_action(x);
}

// prototype declaration
std::string pp(const action_formula& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const action_formula& x)
{
  return out << action_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(action_formula& t1, action_formula& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The value true for action formulas
class true_: public action_formula
{
  public:
    /// \\brief Default constructor X3.
    true_()
      : action_formula(core::detail::default_values::ActTrue)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit true_(const atermpp::aterm& term)
      : action_formula(term)
    {
      assert(core::detail::check_term_ActTrue(*this));
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
  return x.function() == core::detail::function_symbols::ActTrue;
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
  return out << action_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(true_& t1, true_& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The value false for action formulas
class false_: public action_formula
{
  public:
    /// \\brief Default constructor X3.
    false_()
      : action_formula(core::detail::default_values::ActFalse)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit false_(const atermpp::aterm& term)
      : action_formula(term)
    {
      assert(core::detail::check_term_ActFalse(*this));
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
  return x.function() == core::detail::function_symbols::ActFalse;
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
  return out << action_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(false_& t1, false_& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The not operator for action formulas
class not_: public action_formula
{
  public:
    /// \\brief Default constructor X3.
    not_()
      : action_formula(core::detail::default_values::ActNot)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit not_(const atermpp::aterm& term)
      : action_formula(term)
    {
      assert(core::detail::check_term_ActNot(*this));
    }

    /// \\brief Constructor Z14.
    explicit not_(const action_formula& operand)
      : action_formula(atermpp::aterm(core::detail::function_symbol_ActNot(), operand))
    {}

    /// Move semantics
    not_(const not_&) noexcept = default;
    not_(not_&&) noexcept = default;
    not_& operator=(const not_&) noexcept = default;
    not_& operator=(not_&&) noexcept = default;

    const action_formula& operand() const
    {
      return atermpp::down_cast<action_formula>((*this)[0]);
    }
};

/// \\brief Make_not_ constructs a new term into a given address.
/// \\ \param t The reference into which the new not_ is constructed. 
template <class... ARGUMENTS>
inline void make_not_(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_ActNot(), args...);
}

/// \\brief Test for a not expression
/// \\param x A term
/// \\return True if \\a x is a not expression
inline
bool is_not(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::ActNot;
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
  return out << action_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(not_& t1, not_& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The and operator for action formulas
class and_: public action_formula
{
  public:
    /// \\brief Default constructor X3.
    and_()
      : action_formula(core::detail::default_values::ActAnd)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit and_(const atermpp::aterm& term)
      : action_formula(term)
    {
      assert(core::detail::check_term_ActAnd(*this));
    }

    /// \\brief Constructor Z14.
    and_(const action_formula& left, const action_formula& right)
      : action_formula(atermpp::aterm(core::detail::function_symbol_ActAnd(), left, right))
    {}

    /// Move semantics
    and_(const and_&) noexcept = default;
    and_(and_&&) noexcept = default;
    and_& operator=(const and_&) noexcept = default;
    and_& operator=(and_&&) noexcept = default;

    const action_formula& left() const
    {
      return atermpp::down_cast<action_formula>((*this)[0]);
    }

    const action_formula& right() const
    {
      return atermpp::down_cast<action_formula>((*this)[1]);
    }
};

/// \\brief Make_and_ constructs a new term into a given address.
/// \\ \param t The reference into which the new and_ is constructed. 
template <class... ARGUMENTS>
inline void make_and_(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_ActAnd(), args...);
}

/// \\brief Test for a and expression
/// \\param x A term
/// \\return True if \\a x is a and expression
inline
bool is_and(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::ActAnd;
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
  return out << action_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(and_& t1, and_& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The or operator for action formulas
class or_: public action_formula
{
  public:
    /// \\brief Default constructor X3.
    or_()
      : action_formula(core::detail::default_values::ActOr)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit or_(const atermpp::aterm& term)
      : action_formula(term)
    {
      assert(core::detail::check_term_ActOr(*this));
    }

    /// \\brief Constructor Z14.
    or_(const action_formula& left, const action_formula& right)
      : action_formula(atermpp::aterm(core::detail::function_symbol_ActOr(), left, right))
    {}

    /// Move semantics
    or_(const or_&) noexcept = default;
    or_(or_&&) noexcept = default;
    or_& operator=(const or_&) noexcept = default;
    or_& operator=(or_&&) noexcept = default;

    const action_formula& left() const
    {
      return atermpp::down_cast<action_formula>((*this)[0]);
    }

    const action_formula& right() const
    {
      return atermpp::down_cast<action_formula>((*this)[1]);
    }
};

/// \\brief Make_or_ constructs a new term into a given address.
/// \\ \param t The reference into which the new or_ is constructed. 
template <class... ARGUMENTS>
inline void make_or_(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_ActOr(), args...);
}

/// \\brief Test for a or expression
/// \\param x A term
/// \\return True if \\a x is a or expression
inline
bool is_or(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::ActOr;
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
  return out << action_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(or_& t1, or_& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The implication operator for action formulas
class imp: public action_formula
{
  public:
    /// \\brief Default constructor X3.
    imp()
      : action_formula(core::detail::default_values::ActImp)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit imp(const atermpp::aterm& term)
      : action_formula(term)
    {
      assert(core::detail::check_term_ActImp(*this));
    }

    /// \\brief Constructor Z14.
    imp(const action_formula& left, const action_formula& right)
      : action_formula(atermpp::aterm(core::detail::function_symbol_ActImp(), left, right))
    {}

    /// Move semantics
    imp(const imp&) noexcept = default;
    imp(imp&&) noexcept = default;
    imp& operator=(const imp&) noexcept = default;
    imp& operator=(imp&&) noexcept = default;

    const action_formula& left() const
    {
      return atermpp::down_cast<action_formula>((*this)[0]);
    }

    const action_formula& right() const
    {
      return atermpp::down_cast<action_formula>((*this)[1]);
    }
};

/// \\brief Make_imp constructs a new term into a given address.
/// \\ \param t The reference into which the new imp is constructed. 
template <class... ARGUMENTS>
inline void make_imp(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_ActImp(), args...);
}

/// \\brief Test for a imp expression
/// \\param x A term
/// \\return True if \\a x is a imp expression
inline
bool is_imp(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::ActImp;
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
  return out << action_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(imp& t1, imp& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The universal quantification operator for action formulas
class forall: public action_formula
{
  public:
    /// \\brief Default constructor X3.
    forall()
      : action_formula(core::detail::default_values::ActForall)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit forall(const atermpp::aterm& term)
      : action_formula(term)
    {
      assert(core::detail::check_term_ActForall(*this));
    }

    /// \\brief Constructor Z14.
    forall(const data::variable_list& variables, const action_formula& body)
      : action_formula(atermpp::aterm(core::detail::function_symbol_ActForall(), variables, body))
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

    const action_formula& body() const
    {
      return atermpp::down_cast<action_formula>((*this)[1]);
    }
};

/// \\brief Make_forall constructs a new term into a given address.
/// \\ \param t The reference into which the new forall is constructed. 
template <class... ARGUMENTS>
inline void make_forall(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_ActForall(), args...);
}

/// \\brief Test for a forall expression
/// \\param x A term
/// \\return True if \\a x is a forall expression
inline
bool is_forall(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::ActForall;
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
  return out << action_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(forall& t1, forall& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The existential quantification operator for action formulas
class exists: public action_formula
{
  public:
    /// \\brief Default constructor X3.
    exists()
      : action_formula(core::detail::default_values::ActExists)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit exists(const atermpp::aterm& term)
      : action_formula(term)
    {
      assert(core::detail::check_term_ActExists(*this));
    }

    /// \\brief Constructor Z14.
    exists(const data::variable_list& variables, const action_formula& body)
      : action_formula(atermpp::aterm(core::detail::function_symbol_ActExists(), variables, body))
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

    const action_formula& body() const
    {
      return atermpp::down_cast<action_formula>((*this)[1]);
    }
};

/// \\brief Make_exists constructs a new term into a given address.
/// \\ \param t The reference into which the new exists is constructed. 
template <class... ARGUMENTS>
inline void make_exists(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_ActExists(), args...);
}

/// \\brief Test for a exists expression
/// \\param x A term
/// \\return True if \\a x is a exists expression
inline
bool is_exists(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::ActExists;
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
  return out << action_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(exists& t1, exists& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The at operator for action formulas
class at: public action_formula
{
  public:
    /// \\brief Default constructor X3.
    at()
      : action_formula(core::detail::default_values::ActAt)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit at(const atermpp::aterm& term)
      : action_formula(term)
    {
      assert(core::detail::check_term_ActAt(*this));
    }

    /// \\brief Constructor Z14.
    at(const action_formula& operand, const data::data_expression& time_stamp)
      : action_formula(atermpp::aterm(core::detail::function_symbol_ActAt(), operand, time_stamp))
    {}

    /// Move semantics
    at(const at&) noexcept = default;
    at(at&&) noexcept = default;
    at& operator=(const at&) noexcept = default;
    at& operator=(at&&) noexcept = default;

    const action_formula& operand() const
    {
      return atermpp::down_cast<action_formula>((*this)[0]);
    }

    const data::data_expression& time_stamp() const
    {
      return atermpp::down_cast<data::data_expression>((*this)[1]);
    }
};

/// \\brief Make_at constructs a new term into a given address.
/// \\ \param t The reference into which the new at is constructed. 
template <class... ARGUMENTS>
inline void make_at(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_ActAt(), args...);
}

/// \\brief Test for a at expression
/// \\param x A term
/// \\return True if \\a x is a at expression
inline
bool is_at(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::ActAt;
}

// prototype declaration
std::string pp(const at& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const at& x)
{
  return out << action_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(at& t1, at& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The multi action for action formulas
class multi_action: public action_formula
{
  public:
    /// \\brief Default constructor X3.
    multi_action()
      : action_formula(core::detail::default_values::ActMultAct)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit multi_action(const atermpp::aterm& term)
      : action_formula(term)
    {
      assert(core::detail::check_term_ActMultAct(*this));
    }

    /// \\brief Constructor Z14.
    explicit multi_action(const process::action_list& actions)
      : action_formula(atermpp::aterm(core::detail::function_symbol_ActMultAct(), actions))
    {}

    /// Move semantics
    multi_action(const multi_action&) noexcept = default;
    multi_action(multi_action&&) noexcept = default;
    multi_action& operator=(const multi_action&) noexcept = default;
    multi_action& operator=(multi_action&&) noexcept = default;

    const process::action_list& actions() const
    {
      return atermpp::down_cast<process::action_list>((*this)[0]);
    }
};

/// \\brief Make_multi_action constructs a new term into a given address.
/// \\ \param t The reference into which the new multi_action is constructed. 
template <class... ARGUMENTS>
inline void make_multi_action(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_ActMultAct(), args...);
}

/// \\brief Test for a multi_action expression
/// \\param x A term
/// \\return True if \\a x is a multi_action expression
inline
bool is_multi_action(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::ActMultAct;
}

// prototype declaration
std::string pp(const multi_action& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const multi_action& x)
{
  return out << action_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(multi_action& t1, multi_action& t2) noexcept
{
  t1.swap(t2);
}
//--- end generated classes ---//

// template function overloads
std::set<data::variable> find_all_variables(const action_formulas::action_formula& x);

} // namespace mcrl2::action_formulas

#endif // MCRL2_MODAL_FORMULA_ACTION_FORMULA_H
