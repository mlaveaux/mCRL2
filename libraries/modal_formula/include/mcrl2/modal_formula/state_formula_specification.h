// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/modal_formula/state_formula_specification.h
/// \brief add your file description here.

#ifndef MCRL2_MODAL_FORMULA_STATE_FORMULA_SPECIFICATION_H
#define MCRL2_MODAL_FORMULA_STATE_FORMULA_SPECIFICATION_H

#include "mcrl2/modal_formula/state_formula.h"

namespace mcrl2 {

namespace state_formulas {

class state_formula_specification
{
  protected:
    /// \brief The data specification of the specification
    data::data_specification m_data;

    /// \brief The action specification of the specification
    process::action_label_list m_action_labels;

    /// \brief The formula of the specification
    state_formula m_formula;

  public:
    /// \brief Constructor.
    state_formula_specification()
    {}

    /// \brief Constructor of a state formula specification.
    state_formula_specification(const state_formula& formula, const data::data_specification& data = data::data_specification(), const process::action_label_list& action_labels = {})
      : m_data(data),
        m_action_labels(action_labels),
        m_formula(formula)
    {}

    /// \brief Returns the data specification
    /// \return The data specification
    const data::data_specification& data() const
    {
      return m_data;
    }

    /// \brief Returns the data specification
    /// \return The data specification
    data::data_specification& data()
    {
      return m_data;
    }

    /// \brief Returns the action label specification
    /// \return The action label specification
    const process::action_label_list& action_labels() const
    {
      return m_action_labels;
    }

    /// \brief Returns the action label specification
    /// \return The action label specification
    process::action_label_list& action_labels()
    {
      return m_action_labels;
    }

    /// \brief Returns the formula of the state formula specification
    /// \return The formula of the state formula specification
    const state_formula& formula() const
    {
      return m_formula;
    }

    /// \brief Returns the formula of the state formula specification
    /// \return The formula of the state formula specification
    state_formula& formula()
    {
      return m_formula;
    }
};

//--- start generated classes ---//
// prototype declaration
std::string pp(const state_formula_specification& x, bool precedence_aware = true);

/// \\brief Outputs the object to a stream
/// \\param out An output stream
/// \\param x Object x
/// \\return The output stream
inline
std::ostream& operator<<(std::ostream& out, const state_formula_specification& x)
{
  return out << state_formulas::pp(x);
}
//--- end generated classes ---//

} // namespace state_formulas

} // namespace mcrl2

#endif // MCRL2_MODAL_FORMULA_STATE_FORMULA_SPECIFICATION_H
