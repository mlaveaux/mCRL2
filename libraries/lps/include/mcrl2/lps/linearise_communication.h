// Author(s): Jan Friso Groote, Jeroen Keiren
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/lps/linearise_communication.h
/// \brief Apply the rename operator to action summands.

#ifndef MCRL2_LPS_LINEARISE_COMMUNICATION_H
#define MCRL2_LPS_LINEARISE_COMMUNICATION_H

#include "mcrl2/core/detail/print_utility.h"
#include "mcrl2/lps/linearise_allow_block.h"
#include "mcrl2/lps/linearise_utility.h"
#include "mcrl2/lps/stochastic_action_summand.h"
#include "mcrl2/lps/sumelm.h"
#include "mcrl2/lps/detail/configuration.h"
#include "mcrl2/process/process_expression.h"

#include <optional>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/atermpp/aterm_list.cxx"
#else
  import atermpp;
#endif

namespace mcrl2::lps
{

namespace detail
{

/// Get a sorted list of action names that appear in any of the elements of allowlist
/// The result does not contain duplicates.
/// We use this to remove possible communication results that are not relevant because they
/// are guaranteed to be removed from the result after the fact.
inline std::vector<core::identifier_string> get_actions(const process::action_name_multiset_list& allowlist)
{
  std::vector<core::identifier_string> result;
  for (const process::action_name_multiset& l : allowlist)
  {
    const core::identifier_string_list& names = l.names();
    result.insert(result.end(), names.begin(), names.end());
  }
  std::sort(result.begin(), result.end(), action_name_compare());
  const auto it = std::unique(result.begin(), result.end());
  result.erase(it, result.end());

  return result;
}

/// a tuple_list contains pairs of actions (multi-action) and the condition under which this action
/// can occur.
struct tuple_list
{
  std::vector<process::action_list> actions;
  std::vector<data::data_expression> conditions;

  std::size_t size() const
  {
    assert(actions.size() == conditions.size());
    return actions.size();
  }

  tuple_list() = default;

  tuple_list(const std::vector<process::action_list>& actions_, const std::vector<data::data_expression>& conditions_)
      : actions(actions_),
        conditions(conditions_)
  {}

  tuple_list(std::vector<process::action_list>&& actions_, std::vector<data::data_expression>&& conditions_)
      : actions(std::move(actions_)),
        conditions(std::move(conditions_))
  {}
};

/// Extends the list S to S ++ [(a \oplus alpha, c \land c') | (alpha, c') \in T]
///
/// Note that by using move semantics for L, we force the caller to transfer ownership of L to this function,
/// and make it explicit that L should not be used by the caller afterwards.
/// If firstaction == action(), it is not added to the multiactions in L', but the conditions will be strengthened.
/// \pre condition != sort_bool::false_()
inline void addActionCondition(const process::action& a, const data::data_expression& c, tuple_list&& L, tuple_list& S)
{
  assert(c != data::sort_bool::false_()); // It makes no sense to add an action with condition false, as it cannot
                                          // happen anyhow.

  if (a == process::action())
  {
    S.actions.insert(S.actions.end(),
        std::make_move_iterator(L.actions.begin()),
        std::make_move_iterator(L.actions.end()));
  }
  else
  {
    for (process::action_list& m : L.actions)
    {
      m = insert(a, m);
      S.actions.emplace_back(std::move(m));
    }
  }

  if (c == data::sort_bool::true_())
  {
    S.conditions.insert(S.conditions.end(),
        std::make_move_iterator(L.conditions.begin()),
        std::make_move_iterator(L.conditions.end()));
  }
  else
  {
    // Strengthen the conditions in L with condition and append to S.
    for (const data::data_expression& x : L.conditions)
    {
      S.conditions.emplace_back(data::lazy::and_(x, c));
    }
  }
}

/// Data structure to store the communication function more efficiently.
class comm_entry
{
protected:
  /// Left-hand sides of communication expressions
  const std::vector<core::identifier_string_list> m_lhs;

  /// Right-hand sides of communication expressions
  const std::vector<core::identifier_string> m_rhs;

  /// Temporary data using in determining whether communication is allowed.
  /// See usages of the data structure below.
  std::vector<core::identifier_string_list::const_iterator> m_lhs_iters; // offset into lhs
  std::vector<bool> m_match_failed;

  void reset_temporary_data()
  {
    for (std::size_t i = 0; i < size(); ++i)
    {
      m_lhs_iters[i] = m_lhs[i].begin();
      m_match_failed[i] = false;
    }
  }

  /// Check if m is contained in a lhs in the communication entry. (in fact, it must be a prefix).
  /// Returns true if this is the case, false otherwise.
  /// Postcondition: for every i such that m is not contained in lhs[i], match_failed[i] is true.
  /// NB: resets temporary data before performing computations.
  bool match_multiaction(const process::action_list& multi_action)
  {
    assert(std::is_sorted(multi_action.begin(), multi_action.end(), action_compare()));

    reset_temporary_data();

    // m must match a lhs; check every action
    for (const process::action& action : multi_action)
    {
      // check every lhs for actionname
      bool comm_ok = false;
      for (std::size_t i = 0; i < size(); ++i)
      {
        if (m_match_failed[i]) // lhs i does not match
        {
          continue;
        }
        if (m_lhs_iters[i] == m_lhs[i].end()) // lhs cannot match actionname
        {
          m_match_failed[i] = true;
          continue;
        }
        if (action.label().name() != *m_lhs_iters[i])
        {
          // no match
          m_match_failed[i] = true;
        }
        else
        {
          // possible match; on to next action
          ++m_lhs_iters[i];
          comm_ok = true;
        }
      }

      if (!comm_ok) // no (possibly) matching lhs
      {
        return false;
      }
    }

    // There must be an lhs that contains m.
    return true;
  }

  // Initialization of lhs, defined as static function so it can be used in the constructor.
  // Allows lhs to be defined as const.
  static std::vector<core::identifier_string_list> init_lhs(
      const process::communication_expression_list& communications)
  {
    std::vector<core::identifier_string_list> result;
    for (const process::communication_expression& l : communications)
    {
      const core::identifier_string_list& names = l.action_name().names();
      assert(std::is_sorted(names.begin(), names.end(), action_name_compare()));
      result.emplace_back(names.begin(), names.end());
    }
    return result;
  }

  // Initialization of rhs, defined as static function so it can be used in the constructor.
  // Allows rhs to be defined as const.
  static std::vector<core::identifier_string> init_rhs(const process::communication_expression_list& communications)
  {
    std::vector<core::identifier_string> result;
    for (const process::communication_expression& l : communications)
    {
      result.push_back(l.name());
    }
    return result;
  }

public:
  // comm_entries are not copyable.
  comm_entry(const comm_entry&) = delete;
  comm_entry& operator=(const comm_entry&) = delete;

  explicit comm_entry(const process::communication_expression_list& communications)
      : m_lhs(init_lhs(communications)),
        m_rhs(init_rhs(communications)),
        m_lhs_iters(communications.size()),
        m_match_failed(communications.size())
  {}

  ~comm_entry() = default;

  std::size_t size() const
  {
    assert(m_lhs.size() == m_rhs.size());
    assert(m_rhs.size() == m_lhs_iters.size());
    assert(m_lhs_iters.size() == m_match_failed.size());
    return m_lhs.size();
  }

  /// Determine if there exists a communication expression a1|...|an -> b in comm_table
  /// such that m' = a1|...|an , where m' is the multiset of actionnames for multiaction m.
  /// That is, the actions in m consting of actions and data occur in C, such tat a communication
  /// can take place.
  ///
  /// if \exists_{(b,c) \in C} b = \mu(m), return c, otherwise return action_label()
  process::action_label can_communicate(const process::action_list& m)
  {
    assert(std::is_sorted(m.begin(), m.end(), action_compare()));

    process::action_label result; // if no match found, return process::action_label()

    if (match_multiaction(m))
    {
      // There is a lhs that has m as prefix. Find it, and determine if the lhs matches m completely.
      for (std::size_t i = 0; i < size(); ++i)
      {
        // lhs i matches only if comm_table[i] is empty
        if ((!m_match_failed[i]) && m_lhs_iters[i] == m_lhs[i].end())
        {
          if (m_rhs[i] == process::tau())
          {
            throw mcrl2::runtime_error("Cannot linearise a process with a communication operator, containing a "
                                       "communication that results in tau or that has an empty right hand side");
          }
          result = process::action_label(m_rhs[i], m.begin()->label().sorts());
          break;
        }
      }
    }

    return result;
  }

  /// Calculate \exists_{o,c} (\mu(m) \oplus o, c) \in C, with o \subseteq n
  ///
  /// The function calculates whether the actions in m (consisting of actions and data) occur in a left hand side
  /// of a communication a1|...|ak -> b in C (that is, the names of m are a subbag of a1|...|ak), and the actions
  /// a1|...|ak that are not in m are in n. I.e., there is a subbag o of n such that m+o can communicate.
  bool might_communicate(const process::action_list& m, const process::action_list& n)
  {
    assert(std::is_sorted(m.begin(), m.end(), action_compare()));
    assert(std::is_sorted(n.begin(), n.end(), action_compare()));

    /* this function indicates whether the actions in m
       consisting of actions and data occur in C, such that
       a communication might take place (i.e. m is a subbag
       of the lhs of a communication in C).
       if n is not empty, then all actions of a matching communication
       that are not in m should be in n (i.e. there must be a
       subbag o of n such that m+o can communicate. */
    bool result = false;

    if (match_multiaction(m))
    {
      // the rest of actions of lhs that are not in m should be in n
      // rest[i] contains the part of n in which lhs i has to find matching actions
      // if rest[i] cannot match the next remaining action in the left hand side, stored in m_lhs_iters[i], i.e.,
      // rest[i] becomes empty before matching all actions in the lhs, we set it to std::nullopt. N.B. when rest[i]
      // becomes empty after matching all actions in the lhs, rest[i].empty() is a meaningful result: we have a
      // successful match.
      std::vector<std::optional<process::action_list::const_iterator>> rest(size(),
          n.begin()); // pairs of iterator into n; the second element of the pair indicates the end of the range in n.

      // check every lhs
      for (std::size_t i = 0; i < size(); ++i)
      {
        if (m_match_failed[i]) // lhs i did not contain m
        {
          continue;
        }

        // as long as there are still unmatched actions in lhs i...
        while (m_lhs_iters[i] != m_lhs[i].end())
        {
          assert(rest[i] != std::nullopt);
          // .. find them in rest[i]
          if (*rest[i] == n.end()) // no luck
          {
            rest[i] = std::nullopt;
            break;
          }
          // get first action in lhs i
          const core::identifier_string& comm_name = *m_lhs_iters[i];
          core::identifier_string rest_name = (*rest[i])->label().name();
          // find it in rest[i]
          while (comm_name != rest_name)
          {
            ++(*rest[i]);
            if (*rest[i] == n.end()) // no more
            {
              rest[i] = std::nullopt;
              break;
            }
            rest_name = (*rest[i])->label().name();
          }
          if (comm_name != rest_name) // action was not found
          {
            break;
          }

          // action found; try next
          ++(*rest[i]);
          ++m_lhs_iters[i];
        }

        if (rest[i] != std::nullopt) // lhs was found in rest[i]
        {
          result = true;
          break;
        }
      }
    }

    return result;
  }
};

template <typename DataRewriter>
class apply_communication_algorithm
{
public:
  apply_communication_algorithm(const process::action& termination_action,
      DataRewriter& data_rewriter,
      const process::communication_expression_list& communications,
      const process::action_name_multiset_list& allowlist,
      bool is_allow,
      bool is_block)
      : m_terminationAction(termination_action),
        m_data_rewriter(data_rewriter),
        m_communications(sort_communications(communications)),
        m_allowlist(sort_multi_action_labels(allowlist)),
        m_blocked_actions(is_block ? get_actions(allowlist) : std::vector<core::identifier_string>()),
        m_allowed_actions(init_allowed_actions(is_allow, allowlist, termination_action)),
        m_comm_table(m_communications),
        m_is_allow(is_allow),
        m_is_block(is_block)
  {}

  ~apply_communication_algorithm() = default;

  /// Calculate the communication operator applied to a multiaction.
  ///
  /// As the actions in the multiaction can be parameterized with open data expressions, for every subset of
  /// communication expressions that match the multiaction, the result contains a multiaction (in which the particular
  /// subset of communication expressions has been applied), and a corresponding condition that requires that the
  /// parameters are equal for the actions in the left-hand-sides of the matching communication expressions.
  /// Multiactions whose condition is false may be removed.
  ///
  /// This is the function $\overline{\gamma}$ in M. v. Weerdenburg. Calculation of Communication with Open Terms in
  /// GenSpect Process Algebra.
  ///
  /// \param m The multiaction to which the communication operator is applied
  /// \param C The communication expressions to be applied
  /// \param RewriteTerm The rewriter that should be used to simplify the conditions.
  tuple_list apply(const process::action_list& m)
  {
    assert(std::is_sorted(m.begin(), m.end(), action_compare()));
    const process::action_list r;
    return makeMultiActionConditionList_aux(m, r);
  }

  /// Apply the communication composition to a list of action summands.
  void apply(stochastic_action_summand_vector& action_summands,
      deadlock_summand_vector& deadlock_summands,
      bool nosumelm,
      bool nodeltaelimination,
      bool ignore_time)
  {
    assert(!(m_is_allow && m_is_block));

    /* We follow the implementation of Muck van Weerdenburg, described in
       a note: Calculation of communication with open terms. */

    mCRL2log(mcrl2::log::verbose) << (m_is_allow
                                          ? "- calculating the communication operator modulo the allow operator on "
                                      : m_is_block
                                          ? "- calculating the communication operator modulo the block operator on "
                                          : "- calculating the communication operator on ")
                                  << action_summands.size() << " action summands";

    [[maybe_unused]]
    lps_statistics_t lps_statistics_before = get_statistics(action_summands, deadlock_summands);
    // number of summands filtered out after construction of intermediate result (all potential results of communication
    // that may be allowed/are not blocked)
    [[maybe_unused]]
    std::size_t disallowed_summands = 0;      // removed by allow
    [[maybe_unused]]
    std::size_t blocked_summands = 0;         // removed by block
    [[maybe_unused]]
    std::size_t false_condition_summands = 0; // removed because condition is false

    mCRL2log(mcrl2::log::trace) << "Calculating communication operator using a set of " << m_communications.size()
                                << " communication expressions." << std::endl;
    mCRL2log(mcrl2::log::trace) << "Communication expressions: " << std::endl
                                << core::detail::print_set(m_communications) << std::endl;
    mCRL2log(mcrl2::log::trace) << "Allow list: " << std::endl << core::detail::print_set(m_allowlist) << std::endl;

    deadlock_summand_vector resulting_deadlock_summands;
    deadlock_summands.swap(resulting_deadlock_summands);

    const bool inline_allow = m_is_allow || m_is_block;
    if (inline_allow)
    {
      // Inline allow is only supported for ignore_time,
      // for in other cases generation of delta summands cannot be inlined in any simple way.
      assert(!nodeltaelimination && ignore_time);
      deadlock_summands.emplace_back(data::variable_list(), data::sort_bool::true_(), deadlock());
    }

    stochastic_action_summand_vector resulting_action_summands;

    for (const stochastic_action_summand& smmnd : action_summands)
    {
      const data::variable_list& sumvars = smmnd.summation_variables();
      const process::action_list& multiaction = smmnd.multi_action().actions();
      const data::data_expression& time = smmnd.multi_action().time();
      const data::data_expression& condition = smmnd.condition();
      const data::assignment_list& nextstate = smmnd.assignments();
      const stochastic_distribution& dist = smmnd.distribution();

      if (!inline_allow)
      {
        /* Recall a delta summand for every non delta summand.
         * The reason for this is that with communication, the
         * conditions for summands can become much more complex.
         * Many of the actions in these summands are replaced by
         * delta's later on. Due to the more complex conditions it
         * will be hard to remove them. By adding a default delta
         * with a simple condition, makes this job much easier
         * later on, and will in general reduce the number of delta
         * summands in the whole system */

        // Create new list of summand variables containing only those that occur in the condition or the timestamp.
        data::variable_list newsumvars;
        atermpp::make_term_list(
            newsumvars,
            sumvars.begin(),
            sumvars.end(),
            [](const data::variable& v) { return v; },
            [&condition, &time](const data::variable& v)
            { return occursinterm(condition, v) || occursinterm(time, v); });

        resulting_deadlock_summands.emplace_back(newsumvars, condition, deadlock(time));
      }

      /* the multiactionconditionlist is a list containing
         tuples, with a multiaction and the condition,
         expressing whether the multiaction can happen. All
         conditions exclude each other. Furthermore, the list
         is not empty. If no communications can take place,
         the original multiaction is delivered, with condition
         true. */

      // We calculate the communication operator on the multiaction in this single summand. As the actions in the
      // multiaction can be parameterized with open data expressions, for every subset of applicable communication
      // expressions this list in principle contains one summand (unless the condition can be rewritten to false, in
      // which case it is omitted).

      mCRL2log(mcrl2::log::trace) << "Calculating communication on multiaction with " << multiaction.size()
                                  << " actions." << std::endl;
      mCRL2log(mcrl2::log::trace) << "  Multiaction: " << process::pp(multiaction) << std::endl;

      const tuple_list multiactionconditionlist = apply(multiaction);

      mCRL2log(mcrl2::log::trace) << "Calculating communication on multiaction with " << multiaction.size()
                                  << " actions results in " << multiactionconditionlist.size() << " potential summands"
                                  << std::endl;

      for (std::size_t i = 0; i < multiactionconditionlist.size(); ++i)
      {
        const process::action_list& multiaction = multiactionconditionlist.actions[i];

        if (m_is_allow && !allow_(m_allowlist, multiaction, m_terminationAction))
        {
          if constexpr (EnableLineariseStatistics) {
            ++disallowed_summands;
          }

          continue;
        }
        if (m_is_block && encap(m_allowlist, multiaction))
        {
          if constexpr (EnableLineariseStatistics) {
            ++blocked_summands;
          }
          continue;
        }

        const data::data_expression communicationcondition = m_data_rewriter(multiactionconditionlist.conditions[i]);

        const data::data_expression newcondition = m_data_rewriter(data::lazy::and_(condition, communicationcondition));
        stochastic_action_summand new_summand(sumvars,
            newcondition,
            smmnd.multi_action().has_time() ? multi_action(multiaction, smmnd.multi_action().time())
                                            : multi_action(multiaction),
            nextstate,
            dist);
        if (!nosumelm)
        {
          if (sumelm(new_summand))
          {
            new_summand.condition() = m_data_rewriter(new_summand.condition());
          }
        }
        if constexpr (EnableLineariseStatistics)
        {
          if (new_summand.condition() == data::sort_bool::false_())
          {
            ++false_condition_summands;
          }
        }
        
        if (new_summand.condition() != data::sort_bool::false_())
        {
          resulting_action_summands.push_back(new_summand);
        }
      }
    }

    action_summands.swap(resulting_action_summands);

    /* Now the resulting delta summands must be added again */
    if (!inline_allow && !nodeltaelimination)
    {
      for (const deadlock_summand& summand : resulting_deadlock_summands)
      {
        insert_timed_delta_summand(action_summands, deadlock_summands, summand, ignore_time);
      }
    }

    if constexpr (EnableLineariseStatistics)
    {
      lps_statistics_t lps_statistics_after = get_statistics(action_summands, deadlock_summands);
      std::cout << log_comm_application(lps_statistics_before,
          lps_statistics_after,
          disallowed_summands,
          blocked_summands,
          false_condition_summands);
    }

    mCRL2log(mcrl2::log::verbose) << " resulting in " << action_summands.size() << " action summands and "
                                  << deadlock_summands.size() << " delta summands\n";
  }

protected:

  const process::action& m_terminationAction;
  DataRewriter& m_data_rewriter;
  const process::communication_expression_list m_communications;
  const process::action_name_multiset_list m_allowlist;         // This is a list of list of identifierstring.
  const std::vector<core::identifier_string> m_blocked_actions; // used only if m_is_block is set
  const std::vector<core::identifier_string> m_allowed_actions; // used only if m_is_allow is set
  comm_entry m_comm_table;
  const bool
      m_is_allow; // If is_allow or is_block is set, perform inline allow/block filtering. They are mutually exclusive
  const bool m_is_block;

  std::string log_comm_application(const lps_statistics_t& lps_statistics_before,
      const lps_statistics_t& lps_statistics_after,
      const std::size_t disallowed_summands,
      const std::size_t blocked_summands,
      const std::size_t false_condition_summands,
      size_t indent = 0) const
  {
    std::string indent_str(indent, ' ');
    std::ostringstream os;

    os << indent_str << "- operator: comm" << std::endl;

    indent += 2;
    indent_str = std::string(indent, ' ');
    os << indent_str << "number of communication expressions: " << m_communications.size() << std::endl;
    if (m_is_allow)
    {
      os << indent_str << "number of allow expressions: " << m_allowlist.size() << std::endl;
    }
    if (m_is_block)
    {
      os << indent_str << "number of block expressions: " << m_blocked_actions.size() << std::endl;
    }

    os << indent_str << "before:" << std::endl
       << print(lps_statistics_before, indent + 2) << indent_str << "after:" << std::endl
       << print(lps_statistics_after, indent + 2) << indent_str << "removed from intermediate result:" << std::endl;

    indent += 2;
    indent_str = std::string(indent, ' ');
    if (m_is_allow)
    {
      os << indent_str << "disallowed summands: " << disallowed_summands << std::endl;
    }
    if (m_is_block)
    {
      os << indent_str << "blocked summands: " << blocked_summands << std::endl;
    }

    os << indent_str << "summands with false condition: " << false_condition_summands << std::endl;

    return os.str();
  }

  /// Static initialization function to ensure m_allowed_actions can be const.
  static const std::vector<core::identifier_string> init_allowed_actions(bool is_allow,
      const process::action_name_multiset_list& allow_list,
      const process::action& termination_action)
  {
    std::vector<core::identifier_string> result;
    if (is_allow)
    {
      result = get_actions(allow_list);
      // Ensure that the termination action is always allowed, even when it is not an explicit part of
      // the list of allowed actions.
      // We maintains sort order of the vector
      std::vector<core::identifier_string>::iterator it
          = std::upper_bound(result.begin(), result.end(), termination_action.label().name(), action_name_compare());
      result.insert(it, termination_action.label().name());
    }

    return result;
  }

  /// Calculate data expression for pairwise equality of the elements of l1 and l2.
  ///
  /// If the lengths of l1 and l2 differ, or for some index i, the sort l1[i] and l2[i] is different, the pairwise
  /// match is false, otherwise an expression equivalent to \bigwegde_i (l1[i] == l2[i]) is returned.
  data::data_expression pairwise_equal_to(const data::data_expression_list& l1,
      const data::data_expression_list& l2) const
  {
    data::data_expression result = data::sort_bool::true_();

    if (l1.size() != l2.size())
    {
      result = data::sort_bool::false_();
    }

    auto i1 = l1.begin();
    auto i2 = l2.begin();

    while (i1 != l1.end() && i2 != l2.end() && result != data::sort_bool::false_())
    {
      if (i1->sort() != i2->sort())
      {
        result = data::sort_bool::false_();
      }
      else
      {
        result = data::lazy::and_(result, m_data_rewriter(equal_to(*i1++, *i2++)));
      }
    }

    return result;
  }

  /// Determine if the action is allowable, that is, it is part of an allow expression
  /// and it is not part of a block expression
  bool maybe_allowed(const process::action_label& a) const
  {
    assert(std::is_sorted(m_allowed_actions.begin(), m_allowed_actions.end(), action_name_compare()));
    assert(std::is_sorted(m_blocked_actions.begin(), m_blocked_actions.end(), action_name_compare()));
    assert(!m_is_block || m_allowed_actions.empty());
    assert(!m_is_allow || m_blocked_actions.empty());

    return !(m_is_allow || m_is_block)
           || (m_is_allow
               && std::binary_search(m_allowed_actions.begin(),
                   m_allowed_actions.end(),
                   a.name(),
                   action_name_compare()))
           || (m_is_block
               && !std::binary_search(m_blocked_actions.begin(),
                   m_blocked_actions.end(),
                   a.name(),
                   action_name_compare()));
  }

  /// Calculate the communication operator applied to a multiaction.
  ///
  /// This is the function $\overline(\gamma, C, r)$ as described in M. v. Weerdenburg. Calculation of Communication
  /// with Open Terms in GenSpect Process Algebra.
  ///
  /// \param m_first Start of a range of multiactions to which the communication operator should be applied
  /// \param m_last End of a range of multiactions to which the communication operator should be applied
  /// \param C The communication expressions that must be applied to the multiaction
  /// \param r
  /// \param RewriteTerm Data rewriter for simplifying expressions.
  inline tuple_list makeMultiActionConditionList_aux(const process::action_list& m, const process::action_list& r)
  {
    assert(std::is_sorted(m.begin(), m.end(), action_compare()));
    assert(std::is_sorted(r.begin(), r.end(), action_compare()));

    tuple_list S; // result

    // if m = [], then S := { (r, \psi(r, C)) }
    if (m.empty())
    {
      if (r.empty())
      {
        S.conditions.emplace_back(data::sort_bool::true_());
      }
      else
      {
        S.conditions.emplace_back(psi(r));
      }

      // TODO: Why don't we insert r here, as van Weerdenburg writes?
      S.actions.emplace_back();
    }
    else
    {
      // m = [a(d)] \oplus n
      const process::action& a = m.front();
      const process::action_list& m_tail = m.tail();

      // S = \phi(a(d), d, [], n, C, r)
      S = phi({a}, a.arguments(), process::action_list(), m_tail, r);

      // addActionCondition adds a to every multiaction in T; so if a is not part of any allowed action,
      // we can skip this part.
      if (maybe_allowed(a.label()))
      {
        // T = \overline{\gamma}(n, C, [a(d)] \oplus r)
        tuple_list T = makeMultiActionConditionList_aux(m_tail, insert(a, r));

        // S := S \cup \{ (a,true) \oplus t \mid t \in T \}
        // TODO: van Weerdenburg in his note only calculates S := S \cup T. Understand why that is not correct.
        addActionCondition(a, data::sort_bool::true_(), std::move(T), S);
      }
    }
    return S;
  }

  /// Calculate $\phi(m, d, w, n, C, r)$ as described in M. v. Weerdenburg. Calculation of Communication
  /// with Open Terms in GenSpect Process Algebra.
  ///
  /// phi is a function that yields a list of pairs indicating how the actions in m|w|n can communicate.
  /// The pairs contain the resulting multi action and a condition on data indicating when communication can take place.
  /// In the communication all actions of m, none of w and a subset of n can take part in the communication.
  /// d is the data parameter of the communication and comm_table contains a list of multiaction action pairs
  /// indicating possible communications.
  inline tuple_list phi(const process::action_list& m,
      const data::data_expression_list& d,
      const process::action_list& w,
      const process::action_list& n,
      const process::action_list& r)
  {
    assert(std::is_sorted(m.begin(), m.end(), action_compare()));
    assert(std::is_sorted(w.begin(), w.end(), action_compare()));
    assert(std::is_sorted(n.begin(), n.end(), action_compare()));
    assert(std::is_sorted(r.begin(), r.end(), action_compare()));

    tuple_list S;

    // \exists_{o,c} (\mu(m) \oplus o, c) \in C, with o \subseteq n
    if (m_comm_table.might_communicate(m, n))
    {
      if (n.empty()) // b \land n = []
      {
        // There is communication expression whose lhs matches m, result is c.
        const process::action_label c = m_comm_table.can_communicate(m); /* returns process::action_label() if no
                                                                  communication is possible */

        if (c != process::action_label() && maybe_allowed(c))
        {
          // \exists_{(b,c) \in C} b = \mu(m)
          // the result of the communication is part of an allow, not in the block set.
          // Calculate communication for multiaction w.
          // T := \overline{\gamma}(w, C)
          tuple_list T = makeMultiActionConditionList_aux(w, r);

          // S := \{ (c(d) \oplus r, e) \mid (r,e) \in T \}
          addActionCondition(process::action(c, d), data::sort_bool::true_(), std::move(T), S);
        }
      }
      else
      {
        // n = [a(f)] \oplus o
        const process::action& a = n.front();
        const process::action_list& n_tail = n.tail();

        const data::data_expression condition = pairwise_equal_to(d, a.arguments());
        if (condition == data::sort_bool::false_())
        {
          // a(f) cannot take part in communication as the arguments do not match. Move to w and continue with next
          // action
          S = phi(m, d, insert(a, w), n_tail, r);
        }
        else
        {
          tuple_list T = phi(insert(a, m), d, w, n_tail, r);

          S = phi(m, d, insert(a, w), n_tail, r);
          addActionCondition(process::action(), condition, std::move(T), S);
        }
      }
    }

    return S;
  }

  bool xi(const process::action_list& alpha, const process::action_list& beta)
  {
    bool result = false;

    if (beta.empty())
    {
      result = m_comm_table.can_communicate(alpha) != process::action_label();
    }
    else
    {
      const process::action_list alpha_ = insert(beta.front(), alpha);

      if (m_comm_table.can_communicate(alpha_) != process::action_label())
      {
        result = true;
      }
      else
      {
        const process::action_list& beta_tail = beta.tail();
        if (m_comm_table.might_communicate(alpha_, beta_tail))
        {
          result = xi(alpha, beta_tail);
        }

        result = result || xi(alpha, beta_tail);
      }
    }
    return result;
  }

  data::data_expression psi(process::action_list alpha)
  {
    assert(std::is_sorted(alpha.begin(), alpha.end(), action_compare()));
    data::data_expression cond = data::sort_bool::false_();

    process::action_list actl; // used in inner loop.
    while (!alpha.empty())
    {
      process::action_list beta = alpha.tail();

      while (!beta.empty())
      {
        actl = process::action_list();
        actl = insert(alpha.front(), actl);
        actl = insert(beta.front(), actl);
        const process::action_list& beta_tail = beta.tail();
        if (m_comm_table.might_communicate(actl, beta_tail) && xi(actl, beta_tail))
        {
          // sort and remove duplicates??
          cond = data::lazy::or_(cond, pairwise_equal_to(alpha.front().arguments(), beta.front().arguments()));
        }
        beta = beta_tail;
      }

      alpha = alpha.tail();
    }
    cond = data::lazy::not_(cond);
    return cond;
  }
};
} // namespace detail

inline void communicationcomposition(const process::communication_expression_list& communications,
    const process::action_name_multiset_list& allowlist, // This is a list of list of identifierstring.
    const bool is_allow, // If is_allow or is_block is set, perform inline allow/block filtering.
    const bool is_block,
    stochastic_action_summand_vector& action_summands,
    deadlock_summand_vector& deadlock_summands,
    const process::action& terminationAction,
    const bool nosumelm,
    const bool nodeltaelimination,
    const bool ignore_time,
    const std::function<data::data_expression(const data::data_expression&)>& RewriteTerm)

{
  return detail::apply_communication_algorithm(terminationAction,
      RewriteTerm,
      communications,
      allowlist,
      is_allow,
      is_block)
      .apply(action_summands, deadlock_summands, nosumelm, nodeltaelimination, ignore_time);
}

} // namespace mcrl2::lps



#endif // MCRL2_LPS_LINEARISE_COMMUNICATION_H
