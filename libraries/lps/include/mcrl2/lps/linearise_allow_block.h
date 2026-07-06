// Author(s): Jan Friso Groote, Jeroen Keiren
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/lps/linearise_allow_block.h
/// \brief Apply the allow and block operators to summands.

#ifndef MCRL2_LPS_LINEARISE_ALLLOW_BLOCK_H
#define MCRL2_LPS_LINEARISE_ALLLOW_BLOCK_H

#include "mcrl2/atermpp/aterm_list.h"
#include "mcrl2/core/identifier_string.h"
#include "mcrl2/lps/deadlock_summand.h"
#include "mcrl2/lps/detail/configuration.h"
#include "mcrl2/lps/linearise_utility.h"
#include "mcrl2/lps/stochastic_action_summand.h"
#include "mcrl2/process/action_name_multiset.h"
#include "mcrl2/process/process_expression.h"

namespace mcrl2::lps
{

inline std::string log_allow_block_application(const lps_statistics_t& lps_statistics_before,
    const lps_statistics_t& lps_statistics_after,
    const bool is_allow,
    const std::size_t num_allowed_multiactions,
    const std::size_t num_blocked_actions,
    const bool apply_delta_elimination,
    const bool ignore_time,
    size_t indent = 0)
{
  std::string indent_str(indent, ' ');
  std::ostringstream os;

  if (is_allow)
  {
    os << indent_str << "- operator: allow" << std::endl;

    indent += 2;
    indent_str = std::string(indent, ' ');
    os << indent_str << "number of allowed multiactions: " << num_allowed_multiactions << std::endl;
  }
  else
  {
    os << indent_str << "- operator: block" << std::endl;

    indent += 2;
    indent_str = std::string(indent, ' ');
    os << indent_str << "number of blocked actions: " << num_blocked_actions << std::endl;
  }

  os << indent_str << std::boolalpha << "apply delta elimination: " << apply_delta_elimination << std::endl
     << indent_str << "ignore time: " << ignore_time << std::endl
     << indent_str << "before:" << std::endl
     << print(lps_statistics_before, indent + 2) << indent_str << "after:" << std::endl
     << print(lps_statistics_after, indent + 2);

  return os.str();
}

/**************** allow/block *************************************/

namespace detail
{
/// Cache for the allowed actions
/// The cache is used to avoid linear searches in the allow list.
using allow_list_cache = std::unordered_set<core::identifier_string_list>;

/// Calculate the allow_list_cache from the allow_list.
/// This is used to speed up the lookups when determining is a multi-action is allowed.
inline
allow_list_cache make_allow_list_cache(const process::action_name_multiset_list& allowlist)
{
  allow_list_cache result;
  for (const process::action_name_multiset& allow_action : allowlist)
  {
    assert(std::is_sorted(allow_action.names().begin(), allow_action.names().end(), process::action_name_compare()));
    result.insert(allow_action.names());
  }
  return result;
}

/// Calculate the list of action names that appear in a multi-action.
inline
core::identifier_string_list names(const process::action_list& multi_action)
{
  return core::identifier_string_list(
    multi_action.begin(),
    multi_action.end(),
    [](const process::action& a)
    {
      return a.label().name();
    });
}
}

/// \brief Determine if multi_action is allowed by an allow expression in allow_list
///
/// Calculates if the names of the action in multi_action match with an expression in allow_list.
/// If multi_action is the termination action, or multi_action is the empty multiaction, the result is also true.
inline bool allow_(const detail::allow_list_cache& allow_cache,
    const process::action_list& multi_action,
    const process::action& termination_action)
{
  assert(std::is_sorted(multi_action.begin(), multi_action.end()));

  // The empty multiaction and the termination action can never be blocked by allows.
  if (multi_action.empty() || multi_action == process::action_list({termination_action}))
  {
    return true;
  }

  const core::identifier_string_list multi_action_names = detail::names(multi_action);
  return allow_cache.find(multi_action_names) != allow_cache.end();
}

/// \brief Calculate if any of the actions in multiaction is blocked by encap_list.
///
/// \param blocked_actions is a list of action names that are blocked.
/// \param multi_action contains a multiaction a1(d1)|...|an(dn)
/// \returns \f$\exists i: a_i \in encap\_list\f$
template <typename Container>
bool encap(const Container& blocked_actions, const process::action_list& multiaction)
{
  static_assert(std::is_same_v<typename Container::value_type, core::identifier_string>, "Container must contain core::identifier_string");

  assert(std::is_sorted(blocked_actions.begin(), blocked_actions.end(), process::action_name_compare()));
  assert(std::is_sorted(multiaction.begin(), multiaction.end()));

  core::identifier_string_list::const_iterator blocked_actions_it = blocked_actions.begin();
  process::action_list::const_iterator multiaction_it = multiaction.begin();

  while (blocked_actions_it != blocked_actions.end() && multiaction_it != multiaction.end())
  {
    if (*blocked_actions_it == multiaction_it->label().name())
    {
      return true;
    }

    if (process::action_name_compare()(*blocked_actions_it, multiaction_it->label().name()))
    {
      ++blocked_actions_it;
    }
    else
    {
      // The following assertion should hold because blocked_actions_it is not less than or equal to the name
      // of multiaction_it
      assert(process::action_name_compare()(multiaction_it->label().name(), *blocked_actions_it));
      ++multiaction_it;
    }
  }

  return false;
}

/// \brief Calculate if any of the actions in multiaction is blocked by encap_list.
///
/// \param encap_list is a list of action_name_multisets of size 1. Its single element contains all the blocked actions.
/// \param multi_action contains a multiaction a1(d1)|...|an(dn)
/// \returns \f$\exists i: a_i \in encap\_list\f$
inline bool encap(const process::action_name_multiset_list& encaplist, const process::action_list& multiaction)
{
  assert(encaplist.size() == 1);
  assert(std::is_sorted(multiaction.begin(), multiaction.end()));

  const core::identifier_string_list& blocked_actions = encaplist.front().names();
  return encap(blocked_actions, multiaction);
}

/// Calculate the application of the allow or block operator over the action
/// summands.
///
/// \param allowlist1 the allowset. If is_allow is false, the list has length 1,
///                   and its only element contains the blocked actions.
/// \param is_allow determines if we calculate the allow or the block operator
/// \param action_summands the summands over which to compute the operator
/// \param deadlock_summands the deadlock summands tha may be extended.
/// \param termination_action the termination action generated by the linearizer. This action is always allowed.
/// \param ignore_time if true, and nodeltaelimination is false, a true->delta summand is added that subsumes all other
/// deadlock summands \param nodeltaaelimination do not eliminate deadlock summands.
inline void allowblockcomposition(
    const process::action_name_multiset_list& allowlist1, // This is a list of list of identifierstring.
    const bool is_allow,
    stochastic_action_summand_vector& action_summands,
    deadlock_summand_vector& deadlock_summands,
    const process::action& termination_action,
    bool ignore_time,
    bool nodeltaelimination)
{
  // Only keep statistics when these are relevant.
  lps_statistics_t lps_statistics_before = get_statistics(action_summands, deadlock_summands);

  mCRL2log(mcrl2::log::trace) << "Calculating " << ((is_allow) ? "allow" : "block") << " composition using a set of "
                              << ((is_allow) ? allowlist1.size() : allowlist1.front().size())
                              << ((is_allow) ? " allowed multiactions" : " blocked actions") << std::endl;
  mCRL2log(mcrl2::log::trace) << ((is_allow) ? "Allowed multiactions: " : "Blocked actions: ") << std::endl
                              << ((is_allow) ? core::detail::print_set(allowlist1)
                                             : core::detail::print_set(allowlist1.front()))
                              << std::endl;
  /* This function calculates the allow or the block operator,
   depending on whether is_allow is true */

  stochastic_action_summand_vector sourcesumlist;
  action_summands.swap(sourcesumlist);

  deadlock_summand_vector resultdeltasumlist;
  deadlock_summand_vector resultsimpledeltasumlist;
  deadlock_summands.swap(resultdeltasumlist);

  process::action_name_multiset_list allowlist(allowlist1);

  std::size_t sourcesumlist_length = sourcesumlist.size();
  if (sourcesumlist_length > 2 || is_allow) // This condition prevents this message to be printed
  // when performing data elimination. In this case the
  // term delta is linearised, to determine which data
  // is essential for all processes. In these cases a
  // message about the block operator is very confusing.
  {
    mCRL2log(mcrl2::log::verbose) << "- calculating the " << (is_allow ? "allow" : "block") << " operator on "
                                  << sourcesumlist.size() << " action summands and " << resultdeltasumlist.size()
                                  << " delta summands";
  }

  /// Cache the list of allowed actions, to avoid linear searches in the allow list.
  detail::allow_list_cache allow_cache;
  if (is_allow)
  {
    allow_cache = detail::make_allow_list_cache(allowlist);
  }

  /* First add the resulting sums in two separate lists
   one for actions, and one for delta's. The delta's
   are added at the end to the actions, where for
   each delta summand it is determined whether it ought
   to be added, or is superseded by an action or another
   delta summand */
  for (const stochastic_action_summand& smmnd : sourcesumlist)
  {
    const data::variable_list& sumvars = smmnd.summation_variables();
    const process::action_list& multiaction = smmnd.multi_action().actions();
    const data::data_expression& actiontime = smmnd.multi_action().time();
    const data::data_expression& condition = smmnd.condition();

    // Explicitly allow the termination action in any allow.
    if ((is_allow && allow_(allow_cache, multiaction, termination_action))
        || (!is_allow && !encap(allowlist, multiaction)))
    {
      action_summands.push_back(smmnd);
    }
    else if (smmnd.has_time())
    {
      resultdeltasumlist.emplace_back(sumvars, condition, deadlock(actiontime));
    }
    // summand has no time.
    else if (condition == data::sort_bool::true_())
    {
      resultsimpledeltasumlist.emplace_back(sumvars, condition, deadlock());
    }
    else
    {
      resultdeltasumlist.emplace_back(sumvars, condition, deadlock());
    }
  }

  if (nodeltaelimination)
  {
    deadlock_summands.swap(resultsimpledeltasumlist);
    copy(resultdeltasumlist.begin(), resultdeltasumlist.end(), back_inserter(deadlock_summands));
  }
  else if (!ignore_time) /* if a delta summand is added, conditional, timed
                             delta's are subsumed and do not need to be added */
  {
    for (const deadlock_summand& summand : resultsimpledeltasumlist)
    {
      insert_timed_delta_summand(action_summands, deadlock_summands, summand, ignore_time);
    }
    for (const deadlock_summand& summand : resultdeltasumlist)
    {
      insert_timed_delta_summand(action_summands, deadlock_summands, summand, ignore_time);
    }
  }
  else
  {
    // Add a true -> delta
    insert_timed_delta_summand(action_summands,
        deadlock_summands,
        deadlock_summand(data::variable_list(), data::sort_bool::true_(), deadlock()),
        ignore_time);
  }

  if (mCRL2logEnabled(mcrl2::log::verbose) && (sourcesumlist_length > 2 || is_allow))
  {
    mCRL2log(mcrl2::log::verbose) << ", resulting in " << action_summands.size() << " action summands and "
                                  << deadlock_summands.size() << " delta summands\n";
  }

  if constexpr (detail::EnableLineariseStatistics && (sourcesumlist_length > 2 || is_allow))
  {
    // This function is also called when performing data elimination.
    // In that case, there is a block operation that linearizes the term delta.
    // It cannot be correlated to the input process, and the data is confusion.
    lps_statistics_t lps_statistics_after = get_statistics(action_summands, deadlock_summands);

    std::cout << log_allow_block_application(lps_statistics_before,
        lps_statistics_after,
        is_allow,
        (is_allow ? allowlist.size() : 0),
        (is_allow ? 0 : allowlist.front().size()),
        !nodeltaelimination,
        ignore_time);
  }
}

} // namespace mcrl2::lps

#endif // MCRL2_LPS_LINEARISE_ALLLOW_BLOCK_H
