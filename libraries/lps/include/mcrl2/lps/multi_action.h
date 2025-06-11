// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/lps/multi_action.h
/// \brief Multi-action class.

#ifndef MCRL2_LPS_MULTI_ACTION_H
#define MCRL2_LPS_MULTI_ACTION_H

#include "mcrl2/core/print.h"
#include "mcrl2/process/process_expression.h"

namespace mcrl2
{

namespace lps
{

// prototype declaration
bool is_multi_action(const atermpp::aterm& x);

//--- start generated class multi_action ---//
/// \\brief A timed multi-action
class multi_action: public atermpp::aterm
{
  public:


    /// Move semantics
    multi_action(const multi_action&) noexcept = default;
    multi_action(multi_action&&) noexcept = default;
    multi_action& operator=(const multi_action&) noexcept = default;
    multi_action& operator=(multi_action&&) noexcept = default;

    const process::action_list& actions() const
    {
      return atermpp::down_cast<process::action_list>((*this)[0]);
    }

    const data::data_expression& time() const
    {
      return atermpp::down_cast<data::data_expression>((*this)[1]);
    }
//--- start user section multi_action ---//
    /// \brief Constructor
    explicit multi_action(const process::action_list& actions = process::action_list(), 
                          data::data_expression time = data::undefined_real())
      : atermpp::aterm(core::detail::function_symbol_TimedMultAct(), actions, time)
    {
      assert(data::sort_real::is_real(time.sort()));
    }

    /// \brief Constructor.
    /// \param term A term
    explicit multi_action(const atermpp::aterm& term)
      : atermpp::aterm(term)
    {
      assert(core::detail::check_term_TimedMultAct(*this));
    }

    /// \brief Constructor
    explicit multi_action(const process::action& l)
      : multi_action(process::action_list({ l }),data::undefined_real())
    {}

    /// \brief Returns true if time is available.
    /// \return True if time is available.
    bool has_time() const
    {
      return time() != data::undefined_real();
    }

    /// \brief Joins the actions of both multi actions.
    /// \pre The time of both multi actions must be equal.
    multi_action operator+(const multi_action& other) const
    {
      assert(time() == other.time());
      return multi_action(actions() + other.actions(), time());
    }

    /// \brief Returns the multiaction in which the list of actions is sorted. 
    /// \return A multi-action with a sorted list.
    multi_action sort_actions() const
    {
      return multi_action(atermpp::sort_list(actions()),time());
    }

    bool operator==(const multi_action& other) const
    {
      return time()==other.time() && atermpp::sort_list(actions())==atermpp::sort_list(other.actions());
    }

   //--- end user section multi_action ---//
};

/// \\brief Make_multi_action constructs a new term into a given address.
/// \\ \param t The reference into which the new multi_action is constructed. 
template <class... ARGUMENTS>
inline void make_multi_action(atermpp::aterm& t, const ARGUMENTS&... args)
{
  atermpp::make_term_appl(t, core::detail::function_symbol_TimedMultAct(), args...);
}

/// \\brief list of multi_actions
typedef atermpp::term_list<multi_action> multi_action_list;

/// \\brief vector of multi_actions
typedef std::vector<multi_action>    multi_action_vector;

/// \\brief Test for a multi_action expression
/// \\param x A term
/// \\return True if \\a x is a multi_action expression
inline
bool is_multi_action(const atermpp::aterm& x)
{
  return x.function() == core::detail::function_symbols::TimedMultAct;
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
  return out << lps::pp(x);
}

/// \\brief swap overload
inline void swap(multi_action& t1, multi_action& t2)
{
  t1.swap(t2);
}
//--- end generated class multi_action ---//


// template function overloads
lps::multi_action normalize_sorts(const multi_action& x, const data::sort_specification& sortspec);
lps::multi_action translate_user_notation(const lps::multi_action& x);
std::set<data::variable> find_all_variables(const lps::multi_action& x);
std::set<data::variable> find_free_variables(const lps::multi_action& x);

/// \cond INTERNAL_DOCS
namespace detail
{

/// \brief Visits all permutations of the arrays, and calls f for each instance.
/// \pre The range [first, last) contains sorted arrays.
/// \param first Start of a sequence of arrays
/// \param last End of a sequence of arrays
/// \param f A function
template <typename Iter, typename Function>
void forall_permutations(Iter first, Iter last, Function f)
{
  if (first == last)
  {
    f();
    return;
  }
  Iter next = first;
  ++next;
  forall_permutations(next, last, f);
  while (std::next_permutation(first->first, first->second))
  {
    forall_permutations(next, last, f);
  }
}

/// \brief Returns true if the actions in a and b have the same names, and the same sorts.
/// \pre a and b are sorted w.r.t. to the names of the actions.
/// \param a A sequence of actions
/// \param b A sequence of actions
/// \return True if the actions in a and b have the same names, and the same sorts.
inline bool equal_action_signatures(const std::vector<process::action>& a, const std::vector<process::action>& b)
{
  if (a.size() != b.size())
  {
    return false;
  }
  std::vector<process::action>::const_iterator i, j;
  for (i = a.begin(), j = b.begin(); i != a.end(); ++i, ++j)
  {
    if (i->label() != j->label())
    {
      return false;
    }
  }
  return true;
}

/// \brief Compares action labels
struct compare_action_labels
{
  /// \brief Function call operator
  /// \param a An action
  /// \param b An action
  /// \return The function result
  bool operator()(const process::action& a, const process::action& b) const
  {
    return a.label() < b.label();
  }
};

/// \brief Compares action labels and arguments
struct compare_action_label_arguments
{
  /// \brief Function call operator
  /// \param a An action
  /// \param b An action
  /// \return The function result
  bool operator()(const process::action& a, const process::action& b) const
  {
    if (a.label() != b.label())
    {
      return a.label() < b.label();
    }
    return a < b;
  }
};

/// \brief Used for building an expression for the comparison of data parameters.
struct equal_data_parameters_builder
{
  const std::vector<process::action>& a;
  const std::vector<process::action>& b;
  std::set<data::data_expression>& result;

  equal_data_parameters_builder(const std::vector<process::action>& a_,
                                const std::vector<process::action>& b_,
                                std::set<data::data_expression>& result_
                               )
    : a(a_),
      b(b_),
      result(result_)
  {}

  /// \brief Adds the expression 'a == b' to result.
  void operator()()
  {
    std::vector<data::data_expression> v;
    std::vector<process::action>::const_iterator i, j;
    for (i = a.begin(), j = b.begin(); i != a.end(); ++i, ++j)
    {
      data::data_expression_list d1 = i->arguments();
      data::data_expression_list d2 = j->arguments();
      assert(d1.size() == d2.size());
      data::data_expression_list::iterator i1 = d1.begin(), i2 = d2.begin();
      for (     ; i1 != d1.end(); ++i1, ++i2)
      {
        v.push_back(data::lazy::equal_to(*i1, *i2));
      }
    }
    data::data_expression expr = data::lazy::join_and(v.begin(), v.end());
    result.insert(expr);
  }
};

/// \brief Used for building an expression for the comparison of data parameters.
struct not_equal_multi_actions_builder
{
  const std::vector<process::action>& a;
  const std::vector<process::action>& b;
  std::vector<data::data_expression>& result;

  not_equal_multi_actions_builder(const std::vector<process::action>& a_,
                                  const std::vector<process::action>& b_,
                                  std::vector<data::data_expression>& result_
                                 )
    : a(a_),
      b(b_),
      result(result_)
  {}

  /// \brief Adds the expression 'a == b' to result.
  void operator()()
  {
    using namespace data::lazy;

    std::vector<data::data_expression> v;
    std::vector<process::action>::const_iterator i, j;
    for (i = a.begin(), j = b.begin(); i != a.end(); ++i, ++j)
    {
      data::data_expression_list d1 = i->arguments();
      data::data_expression_list d2 = j->arguments();
      assert(d1.size() == d2.size());
      data::data_expression_list::iterator i1=d1.begin(), i2=d2.begin();
      for (   ; i1 != d1.end(); ++i1, ++i2)
      {
        v.push_back(data::not_equal_to(*i1, *i2));
      }
    }
    result.push_back(data::lazy::join_or(v.begin(), v.end()));
  }
};

} // namespace detail
/// \endcond

/// \brief Returns a data expression that expresses under which conditions the
/// multi actions a and b are equal. The multi actions may contain free variables.
/// \param a A sequence of actions
/// \param b A sequence of actions
/// \return Necessary conditions for the equality of a and b
inline data::data_expression equal_multi_actions(const multi_action& a, const multi_action& b)
{
#ifdef MCRL2_EQUAL_MULTI_ACTIONS_DEBUG
  mCRL2log(debug) << "\n<equal multi actions>" << std::endl;
  mCRL2log(debug) << "a = " << process::pp(a.actions()) << std::endl;
  mCRL2log(debug) << "b = " << process::pp(b.actions()) << std::endl;
#endif
  using namespace data::lazy;

  // make copies of a and b and sort them
  std::vector<process::action> va(a.actions().begin(), a.actions().end()); // protection not needed
  std::vector<process::action> vb(b.actions().begin(), b.actions().end()); // protection not needed
  std::sort(va.begin(), va.end(), detail::compare_action_label_arguments());
  std::sort(vb.begin(), vb.end(), detail::compare_action_label_arguments());

  if (!detail::equal_action_signatures(va, vb))
  {
#ifdef MCRL2_EQUAL_MULTI_ACTIONS_DEBUG
    mCRL2log(debug) << "different action signatures detected!" << std::endl;
    mCRL2log(debug) << "a = " << process::action_list(va.begin(), va.end()) << std::endl;
    mCRL2log(debug) << "b = " << process::action_list(vb.begin(), vb.end()) << std::endl;
#endif
    return data::sort_bool::false_();
  }

  // compute the intervals of a with equal names
  typedef std::vector<process::action>::iterator action_iterator;
  std::vector<std::pair<action_iterator, action_iterator> > intervals;
  action_iterator first = va.begin();
  while (first != va.end())
  {
    action_iterator next = std::upper_bound(first, va.end(), *first, detail::compare_action_labels());
    intervals.push_back(std::make_pair(first, next));
    first = next;
  }

  std::set<data::data_expression> z;
  detail::equal_data_parameters_builder f(va, vb, z);
  detail::forall_permutations(intervals.begin(), intervals.end(), f);
  data::data_expression result = data::lazy::join_or(z.begin(), z.end());
  return result;
}

/// \brief Returns a pbes expression that expresses under which conditions the
/// multi actions a and b are not equal. The multi actions may contain free variables.
/// \param a A sequence of actions
/// \param b A sequence of actions
/// \return Necessary conditions for the inequality of a and b
inline data::data_expression not_equal_multi_actions(const multi_action& a, const multi_action& b)
{
  using namespace data::lazy;

  // make copies of a and b and sort them
  std::vector<process::action> va(a.actions().begin(), a.actions().end());
  std::vector<process::action> vb(b.actions().begin(), b.actions().end());
  std::sort(va.begin(), va.end(), detail::compare_action_label_arguments());
  std::sort(vb.begin(), vb.end(), detail::compare_action_label_arguments());

  if (!detail::equal_action_signatures(va, vb))
  {
    return data::sort_bool::true_();
  }

  // compute the intervals of a with equal names
  typedef std::vector<process::action>::iterator action_iterator;
  std::vector<std::pair<action_iterator, action_iterator> > intervals;
  action_iterator first = va.begin();
  while (first != va.end())
  {
    action_iterator next = std::upper_bound(first, va.end(), *first, detail::compare_action_labels());
    intervals.push_back(std::make_pair(first, next));
    first = next;
  }
  std::vector<data::data_expression> z;
  detail::not_equal_multi_actions_builder f(va, vb, z);
  detail::forall_permutations(intervals.begin(), intervals.end(), f);
  data::data_expression result = data::lazy::join_and(z.begin(), z.end());
  return result;
}

} // namespace lps

} // namespace mcrl2

namespace std
{
/// \brief specialization of the standard std::hash function for an action_label_string.
template<>
struct hash< mcrl2::lps::multi_action >
{
  std::size_t operator()(const mcrl2::lps::multi_action& ma) const
  {
    std::hash<atermpp::aterm> hasher;
    return hasher(ma.actions()) ^ (hasher(ma.time())<<1);
  }
};

} // namespace std

#endif // MCRL2_LPS_MULTI_ACTION_H
