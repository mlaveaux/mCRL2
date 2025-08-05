// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/modal_formula/regular_formula.h
/// \brief Add your file description here.

#ifndef MCRL2_MODAL_FORMULA_REGULAR_FORMULA_H
#define MCRL2_MODAL_FORMULA_REGULAR_FORMULA_H

#include "mcrl2/modal_formula/action_formula.h"

namespace mcrl2::regular_formulas
{

//--- start generated classes ---//
/// \\brief A regular formula
class regular_formula: public atermpp::aterm
{
  public:
    /// \\brief Default constructor X3.
    regular_formula()
      : atermpp::aterm(core::detail::default_values::RegFrm)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit regular_formula(const atermpp::aterm& term)
      : atermpp::aterm(term)
    {
      assert(core::detail::check_rule_RegFrm(*this));
    }

    /// \\brief Constructor Z6.
    regular_formula(const action_formulas::action_formula& x)
      : atermpp::aterm(x)
    {}

    /// \\brief Constructor Z6.
    regular_formula(const data::data_expression& x)
      : atermpp::aterm(x)
    {}

    /// Move semantics
    regular_formula(const regular_formula&) noexcept = default;
    regular_formula(regular_formula&&) noexcept = default;
    regular_formula& operator=(const regular_formula&) noexcept = default;
    regular_formula& operator=(regular_formula&&) noexcept = default;
};

/// \\brief list of regular_formulas
using regular_formula_list = atermpp::term_list<regular_formula>;

/// \\brief vector of regular_formulas
using regular_formula_vector = std::vector<regular_formula>;

// prototypes
inline bool is_seq(const atermpp::aterm& x);
inline bool is_alt(const atermpp::aterm& x);
inline bool is_trans(const atermpp::aterm& x);
inline bool is_trans_or_nil(const atermpp::aterm& x);
inline bool is_untyped_regular_formula(const atermpp::aterm& x);

/// \\brief Test for a regular_formula expression
/// \\param x A term
/// \\return True if \\a x is a regular_formula expression
inline
bool is_regular_formula(const atermpp::aterm& x)
{
  return data::is_data_expression(x) ||
         action_formulas::is_action_formula(x) ||
         regular_formulas::is_seq(x) ||
         regular_formulas::is_alt(x) ||
         regular_formulas::is_trans(x) ||
         regular_formulas::is_trans_or_nil(x) ||
         regular_formulas::is_untyped_regular_formula(x);
}

// prototype declaration
std::string pp(const regular_formula& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const regular_formula& x)
{
  return out << regular_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(regular_formula& t1, regular_formula& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The seq operator for regular formulas
class seq: public regular_formula
{
  public:
    /// \\brief Default constructor X3.
    seq()
      : regular_formula(core::detail::default_values::RegSeq)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit seq(const atermpp::aterm& term)
      : regular_formula(term)
    {
      assert(core::detail::check_term_RegSeq(*this));
    }

    /// \\brief Constructor Z14.
    seq(const regular_formula& left, const regular_formula& right)
      : regular_formula(atermpp::aterm(core::detail::function_symbol_RegSeq(), left, right))
    {}

    /// Move semantics
    seq(const seq&) noexcept = default;
    seq(seq&&) noexcept = default;
    seq& operator=(const seq&) noexcept = default;
    seq& operator=(seq&&) noexcept = default;

    const regular_formula& left() const
    {
      return atermpp::down_cast<regular_formula>((*this)[0]);
    }

    const regular_formula& right() const
    {
      return atermpp::down_cast<regular_formula>((*this)[1]);
    }
};

/// \\brief Make_seq constructs a new term into a given address.
/// \\ \param t The reference into which the new seq is constructed. 
template <class... ARGUMENTS>
inline void make_seq(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_RegSeq(), args...);
}

/// \\brief Test for a seq expression
/// \\param x A term
/// \\return True if \\a x is a seq expression
inline
bool is_seq(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::RegSeq;
}

// prototype declaration
std::string pp(const seq& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const seq& x)
{
  return out << regular_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(seq& t1, seq& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The alt operator for regular formulas
class alt: public regular_formula
{
  public:
    /// \\brief Default constructor X3.
    alt()
      : regular_formula(core::detail::default_values::RegAlt)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit alt(const atermpp::aterm& term)
      : regular_formula(term)
    {
      assert(core::detail::check_term_RegAlt(*this));
    }

    /// \\brief Constructor Z14.
    alt(const regular_formula& left, const regular_formula& right)
      : regular_formula(atermpp::aterm(core::detail::function_symbol_RegAlt(), left, right))
    {}

    /// Move semantics
    alt(const alt&) noexcept = default;
    alt(alt&&) noexcept = default;
    alt& operator=(const alt&) noexcept = default;
    alt& operator=(alt&&) noexcept = default;

    const regular_formula& left() const
    {
      return atermpp::down_cast<regular_formula>((*this)[0]);
    }

    const regular_formula& right() const
    {
      return atermpp::down_cast<regular_formula>((*this)[1]);
    }
};

/// \\brief Make_alt constructs a new term into a given address.
/// \\ \param t The reference into which the new alt is constructed. 
template <class... ARGUMENTS>
inline void make_alt(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_RegAlt(), args...);
}

/// \\brief Test for a alt expression
/// \\param x A term
/// \\return True if \\a x is a alt expression
inline
bool is_alt(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::RegAlt;
}

// prototype declaration
std::string pp(const alt& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const alt& x)
{
  return out << regular_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(alt& t1, alt& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The trans operator for regular formulas
class trans: public regular_formula
{
  public:
    /// \\brief Default constructor X3.
    trans()
      : regular_formula(core::detail::default_values::RegTrans)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit trans(const atermpp::aterm& term)
      : regular_formula(term)
    {
      assert(core::detail::check_term_RegTrans(*this));
    }

    /// \\brief Constructor Z14.
    explicit trans(const regular_formula& operand)
      : regular_formula(atermpp::aterm(core::detail::function_symbol_RegTrans(), operand))
    {}

    /// Move semantics
    trans(const trans&) noexcept = default;
    trans(trans&&) noexcept = default;
    trans& operator=(const trans&) noexcept = default;
    trans& operator=(trans&&) noexcept = default;

    const regular_formula& operand() const
    {
      return atermpp::down_cast<regular_formula>((*this)[0]);
    }
};

/// \\brief Make_trans constructs a new term into a given address.
/// \\ \param t The reference into which the new trans is constructed. 
template <class... ARGUMENTS>
inline void make_trans(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_RegTrans(), args...);
}

/// \\brief Test for a trans expression
/// \\param x A term
/// \\return True if \\a x is a trans expression
inline
bool is_trans(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::RegTrans;
}

// prototype declaration
std::string pp(const trans& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const trans& x)
{
  return out << regular_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(trans& t1, trans& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief The 'trans or nil' operator for regular formulas
class trans_or_nil: public regular_formula
{
  public:
    /// \\brief Default constructor X3.
    trans_or_nil()
      : regular_formula(core::detail::default_values::RegTransOrNil)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit trans_or_nil(const atermpp::aterm& term)
      : regular_formula(term)
    {
      assert(core::detail::check_term_RegTransOrNil(*this));
    }

    /// \\brief Constructor Z14.
    explicit trans_or_nil(const regular_formula& operand)
      : regular_formula(atermpp::aterm(core::detail::function_symbol_RegTransOrNil(), operand))
    {}

    /// Move semantics
    trans_or_nil(const trans_or_nil&) noexcept = default;
    trans_or_nil(trans_or_nil&&) noexcept = default;
    trans_or_nil& operator=(const trans_or_nil&) noexcept = default;
    trans_or_nil& operator=(trans_or_nil&&) noexcept = default;

    const regular_formula& operand() const
    {
      return atermpp::down_cast<regular_formula>((*this)[0]);
    }
};

/// \\brief Make_trans_or_nil constructs a new term into a given address.
/// \\ \param t The reference into which the new trans_or_nil is constructed. 
template <class... ARGUMENTS>
inline void make_trans_or_nil(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_RegTransOrNil(), args...);
}

/// \\brief Test for a trans_or_nil expression
/// \\param x A term
/// \\return True if \\a x is a trans_or_nil expression
inline
bool is_trans_or_nil(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::RegTransOrNil;
}

// prototype declaration
std::string pp(const trans_or_nil& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const trans_or_nil& x)
{
  return out << regular_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(trans_or_nil& t1, trans_or_nil& t2) noexcept
{
  t1.swap(t2);
}


/// \\brief An untyped regular formula or action formula
class untyped_regular_formula: public regular_formula
{
  public:
    /// \\brief Default constructor X3.
    untyped_regular_formula()
      : regular_formula(core::detail::default_values::UntypedRegFrm)
    {}

    /// \\brief Constructor Z9.
    /// \\param term A term
    explicit untyped_regular_formula(const atermpp::aterm& term)
      : regular_formula(term)
    {
      assert(core::detail::check_term_UntypedRegFrm(*this));
    }

    /// \\brief Constructor Z14.
    untyped_regular_formula(const core::identifier_string& name, const regular_formula& left, const regular_formula& right)
      : regular_formula(atermpp::aterm(core::detail::function_symbol_UntypedRegFrm(), name, left, right))
    {}

    /// \\brief Constructor Z2.
    untyped_regular_formula(const std::string& name, const regular_formula& left, const regular_formula& right)
      : regular_formula(atermpp::aterm(core::detail::function_symbol_UntypedRegFrm(), core::identifier_string(name), left, right))
    {}

    /// Move semantics
    untyped_regular_formula(const untyped_regular_formula&) noexcept = default;
    untyped_regular_formula(untyped_regular_formula&&) noexcept = default;
    untyped_regular_formula& operator=(const untyped_regular_formula&) noexcept = default;
    untyped_regular_formula& operator=(untyped_regular_formula&&) noexcept = default;

    const core::identifier_string& name() const
    {
      return atermpp::down_cast<core::identifier_string>((*this)[0]);
    }

    const regular_formula& left() const
    {
      return atermpp::down_cast<regular_formula>((*this)[1]);
    }

    const regular_formula& right() const
    {
      return atermpp::down_cast<regular_formula>((*this)[2]);
    }
};

/// \\brief Make_untyped_regular_formula constructs a new term into a given address.
/// \\ \param t The reference into which the new untyped_regular_formula is constructed. 
template <class... ARGUMENTS>
inline void make_untyped_regular_formula(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_UntypedRegFrm(), args...);
}

/// \\brief Test for a untyped_regular_formula expression
/// \\param x A term
/// \\return True if \\a x is a untyped_regular_formula expression
inline
bool is_untyped_regular_formula(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::UntypedRegFrm;
}

// prototype declaration
std::string pp(const untyped_regular_formula& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const untyped_regular_formula& x)
{
  return out << regular_formulas::pp(x);
}

/// \\brief swap overload
inline void swap(untyped_regular_formula& t1, untyped_regular_formula& t2) noexcept
{
  t1.swap(t2);
}
//--- end generated classes ---//

} // namespace mcrl2::regular_formulas



#endif // MCRL2_MODAL_FORMULA_REGULAR_FORMULA_H
