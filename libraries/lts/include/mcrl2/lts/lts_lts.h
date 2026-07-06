// Author(s): Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

/** \file lts_lts.h
 *
 * \brief This file contains a class that contains labelled transition systems in lts (mcrl2) format.
 * \details A labelled transition system in lts/mcrl2 format is a transition system
 * with as state labels vectors of strings, and as transition labels strings.
 * \author Jan Friso Groote
 */


#ifndef MCRL2_LTS_LTS_MCRL2_H
#define MCRL2_LTS_LTS_MCRL2_H

#include "mcrl2/lps/parse.h"
#include "mcrl2/lps/print.h"
#include "mcrl2/lps/state.h"
#include "mcrl2/lps/probabilistic_data_expression.h"
#include "mcrl2/lts/probabilistic_lts.h"

namespace mcrl2::lts
{

/** \brief This class contains state labels for an labelled transition system in .lts format.
    \details A state label in .lts format consists of lists of balanced tree of data expressions.
             These represent sets of state vectors. The reason for the sets is that states can
             be merged by operations on state spaces, and if so, the sets of labels can easily
             be joined.
*/
class state_label_lts : public atermpp::term_list< lps::state >
{
  public:
    using super = atermpp::term_list<lps::state>;

    /** \brief Default constructor.
    */
    state_label_lts() = default;

    /** \brief Copy constructor. */
    state_label_lts(const state_label_lts& )=default;

    /** \brief Copy assignment. */
    state_label_lts& operator=(const state_label_lts& )=default;

    /** \brief Construct a single state label out of the elements in a container.
    */
    template < class CONTAINER >
    explicit state_label_lts(const CONTAINER& l)
    {
      static_assert(std::is_same_v<typename CONTAINER::value_type, data::data_expression>,
          "Value type must be a data_expression");
      this->push_front(lps::state(l.begin(),l.size()));
    }

    /** \brief Construct a state label out of a balanced tree of data expressions, representing a state label.
    */
    explicit state_label_lts(const lps::state& l)
    {
      this->push_front(l);
    }

    /** \brief Construct a state label out of list of balanced trees of data expressions, representing a state label.
    */
    explicit state_label_lts(const super& l)
      : super(l)
    {
    }

    /** \brief An operator to concatenate two state labels.
        \details Is optimal whenever |l| is smaller than the left operand, i.e. |l| < |this|.
     */
    state_label_lts operator+(const state_label_lts& l) const
    {
      // The order of the state labels should not matter. For efficiency the elements of l are inserted in front of the aterm_list.
      state_label_lts result(*this);

      for (const lps::state& el : l)
      {
        result.push_front(el);
      }

      return result;
    }

    /** \brief Create a state label consisting of a number as the only list element.
    */
    static state_label_lts number_to_label(const std::size_t n)
    {
      return state_label_lts(lps::state(mcrl2::data::sort_nat::nat(n)));
    }
};

/** \brief Pretty print a state value of this LTS.
    \details The label is printed as (t1,...,tn) if it consists of a single label.
             Otherwise the labels are printed with square brackets surrounding it.
    \param[in] label  The state value to pretty print.
    \return           The pretty-printed representation of value. */

inline std::string pp(const state_label_lts& label)
{
  std::string s;
  if (label.size()!=1)
  {
    s += "[";
  }
  bool first=true;
  for(const lps::state& l: label)
  {
    if (!first)
    {
      s += ", ";
    }
    first = false;
    s += "(";
    for (lps::state::iterator i=l.begin(); i!=l.end(); ++i)
    {
      if (i!=l.begin())
      {
        s += ",";
      }
      s += data::pp(*i);
    }
    s += ")";
  }
  if (label.size()!=1)
  {
    s += "]";
  }
  return s;
}

/** \brief A class containing the values for action labels for the .lts format.
    \details An action label is a multi_action. */
class action_label_lts: public mcrl2::lps::multi_action
{
  public:

    /** \brief Default constructor. */
    action_label_lts() = default;

    /** \brief Copy constructor. */
    action_label_lts(const action_label_lts& )=default;

    /** \brief Copy assignment. */
    action_label_lts& operator=(const action_label_lts& )=default;

    /** \brief Constructor. */
    explicit action_label_lts(const mcrl2::lps::multi_action& a)
     : mcrl2::lps::multi_action(a)
    {
    }

    /** \brief Hide the actions with labels in tau_actions.
    */
    void hide_actions(const std::vector<std::string>& tau_actions)
    {
      process::action_vector new_multi_action;
      for (const process::action& a: this->actions()) 
      {
        if (std::find(tau_actions.begin(),tau_actions.end(),
                      std::string(a.label().name()))==tau_actions.end())  // this action must not be hidden.
        {
          new_multi_action.push_back(a);
        }
      }
      *this= action_label_lts(lps::multi_action(process::action_list(new_multi_action.begin(), new_multi_action.end()),time()));
    }

    /* \brief The action label that represents the internal action.
    */
    static const action_label_lts& tau_action()
    {
      static action_label_lts tau_action=action_label_lts();
      return tau_action;
    }

};

/** \brief Print the action label to string. */
inline std::string pp(const action_label_lts& l)
{
  return lps::pp(mcrl2::lps::multi_action(l));
}

/** \brief Parse a string into an action label.
    \details The string is typechecked against the data specification and
             list of declared actions. If parsing or type checking fails, an
             mcrl2::runtime_error is thrown.
    \param[in] multi_action_string The string to be parsed.
    \param[in] data_spec The data specification used for parsing.
    \param[in] typechecker The type checker used for parsing.
    \return The parsed and type checked multi action. */
inline action_label_lts parse_lts_action(
  const std::string& multi_action_string,
  const data::data_specification& data_spec,
  lps::multi_action_type_checker& typechecker)
{
  std::string l(multi_action_string);
  lps::multi_action al;
  // Find an @ symbol in the action, in which case it is timed.
  size_t position_of_at=l.find_first_of('@');
  std::string time_expression_string;
  if (position_of_at!=std::string::npos)
  {
    // The at symbol is found. It is a timed action.
    time_expression_string=l.substr(position_of_at+1,l.size());
    l=l.substr(0,position_of_at-1);
  }
  try
  {
    al=lps::parse_multi_action(l, typechecker, data_spec);
  }
  catch (mcrl2::runtime_error& e)
  {
    throw mcrl2::runtime_error("Parse error in action label " + multi_action_string + ".\n" + e.what());
  }

  if (time_expression_string.size()>0)
  {
    try
    {
      data::data_expression time=parse_data_expression(time_expression_string,data::variable_list(),data_spec);
      // Translate the sort of time to a real.
      if (time.sort()==data::sort_pos::pos())
      {
#ifdef MCRL2_ENABLE_MACHINENUMBERS
        time = data::sort_nat::pos2nat(time);
#else
        time = data::sort_nat::cnat(time);
#endif
      }
      if (time.sort()==data::sort_nat::nat())
      {
        time = data::sort_int::cint(time);
      }
      if (time.sort()==data::sort_int::int_())
      {
        time = data::sort_real::creal(time, data::sort_pos::c1());
      }
      if (time.sort()!=data::sort_real::real_())
      {
        throw mcrl2::runtime_error("The time is not of sort Pos, Nat, Int or Real\n");
      }
      return action_label_lts(lps::multi_action(al.actions(),time));
    }
    catch (mcrl2::runtime_error& e)
    {
      throw mcrl2::runtime_error("Parse error in the time expression " + multi_action_string + ".\n" + e.what());
    }
  }
  return action_label_lts(al);
}

namespace detail
{

/** \brief a base class for lts_lts_t and probabilistic_lts_t.
 */
class lts_lts_base
{
  protected:
    data::data_specification m_data_spec;
    data::variable_list m_parameters;
    process::action_label_list m_action_decls;

  public:
    /// \brief Default constructor
    lts_lts_base() = default;

    /// \brief Standard equality function;
    bool operator==(const lts_lts_base& other) const
    {
      return m_data_spec==other.m_data_spec &&
             m_parameters==other.m_parameters &&
             m_action_decls==other.m_action_decls;
    }

    void swap(lts_lts_base& l) noexcept
    {
      const data::data_specification auxd=m_data_spec;
      m_data_spec=l.m_data_spec;
      l.m_data_spec=auxd;
      m_parameters.swap(l.m_parameters);
      m_action_decls.swap(l.m_action_decls);
    }

    /** \brief Yields the type of this lts, in this case lts_lts. */
    static lts_type type()
    {
      return lts_lts;
    }

    /** \brief Returns the mCRL2 data specification of this LTS.
    */
    const data::data_specification& data() const
    {
      return m_data_spec;
    }

    /** \brief Return action label declarations stored in this LTS.
    */
    const process::action_label_list& action_label_declarations() const
    {
      return m_action_decls;
    }

    /** \brief Set the action label information for this LTS.
     * \param[in] decls  The action labels to be set in this lts.
    */
    void set_action_label_declarations(const process::action_label_list& decls)
    {
      m_action_decls=decls;
    }

    /** \brief Set the mCRL2 data specification of this LTS.
     * \param[in] spec  The mCRL2 data specification for this LTS.
    */
    void set_data(const data::data_specification& spec)
    {
      m_data_spec=spec;
    }

    /** \brief Return the process parameters stored in this LTS.
    */
    const data::variable_list& process_parameters() const
    {
      return m_parameters;
    }

    /** \brief Returns the i-th parameter of the state vectors stored in this LTS.
     * \return The state parameters stored in this LTS.
    */
    const data::variable& process_parameter(std::size_t i) const
    {
      assert(i<m_parameters.size());
      data::variable_list::const_iterator j=m_parameters.begin();
      for(std::size_t c=0; c!=i; ++j, ++c)
      {}
      return *j;
    }


    /** \brief Set the state parameters for this LTS.
     * \param[in] params  The state parameters for this lts.
    */
    void set_process_parameters(const data::variable_list& params)
    {
      m_parameters=params;
    }

};

} // namespace detail


/** \brief This class contains labelled transition systems in .lts format.
    \details In this .lts format, an action label is a multi action, and a
           state label is an expression of the form STATE(t1,...,tn) where
           ti are data expressions.
*/
class lts_lts_t : public lts< state_label_lts, action_label_lts, detail::lts_lts_base >
{
  public:
    /** \brief Creates an object containing no information. */
    lts_lts_t() = default;

    /** \brief Load the labelled transition system from file.
     *  \details If the filename is empty, the result is read from stdout.
     *  \param[in] filename Name of the file to which this lts is written.
     */
    void load(const std::string& filename);

    /** \brief Save the labelled transition system to file.
     *  \details If the filename is empty, the result is read from stdin.
     *  \param[in] filename Name of the file from which this lts is read.
     */
    void save(const std::string& filename) const;
};

/** \brief This class contains probabilistic labelled transition systems in .lts format.
    \details In this .lts format, an action label is a multi action, and a
           state label is an expression of the form STATE(t1,...,tn) where
           ti are data expressions.
*/
class probabilistic_lts_lts_t :
         public probabilistic_lts< state_label_lts,
                                   action_label_lts,
                                   probabilistic_state<std::size_t, lps::probabilistic_data_expression>,
                                   detail::lts_lts_base >
{
  public:
    /** \brief Creates an object containing no information. */
    probabilistic_lts_lts_t() = default;

    /** \brief Load the labelled transition system from file.
     *  \details If the filename is empty, the result is read from stdout.
     *  \param[in] filename Name of the file to which this lts is written.
     */
    void load(const std::string& filename);

    /** \brief Save the labelled transition system to file.
     *  \details If the filename is empty, the result is read from stdin.
     *  \param[in] filename Name of the file from which this lts is read.
     */
    void save(const std::string& filename) const;
};
} // namespace mcrl2::lts

namespace std
{

/// \brief specialization of the standard std::hash function for an action_label_string.
template<>
struct hash< mcrl2::lts::action_label_lts >
{
  std::size_t operator()(const mcrl2::lts::action_label_lts& as) const
  {
    hash<mcrl2::lps::multi_action> hasher;
    return hasher(as);
  }
};

} // namespace std


#endif
