// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/lts/state_space_generator.h
/// \brief add your file description here.

#ifndef MCRL2_LTS_STATE_SPACE_GENERATOR_H
#define MCRL2_LTS_STATE_SPACE_GENERATOR_H

#include <forward_list>

#include "mcrl2/lps/explorer.h"
#include "mcrl2/lts/trace.h"

namespace mcrl2::lts 
{

inline
std::ostream& operator<<(std::ostream& out, const lps::state& s)
{
  return out << atermpp::pp(s);
}

inline
const lps::state& first_state(const lps::state& s)
{
  return s;
}

inline
const lps::state& first_state(const lps::stochastic_state& s)
{
  return s.states.front();
}

namespace detail 
{

inline
bool save_trace(
  class trace& tr,
  const std::string& filename
)
{
  try
  {
    tr.save(filename);
    mCRL2log(log::info) << " and saved trace to '" << filename << "'";
    return true;
  }
  catch(...)
  {
    mCRL2log(log::info) << ", but its trace could not be saved to '" << filename << "'";
  }
  return false;
}

inline
void save_traces(
  class trace& tr,
  const std::string& filename1,
  class trace& tr2,
  const std::string& filename2
)
{
  try
  {
    tr.save(filename1);
    tr2.save(filename2);
    mCRL2log(log::info) << " and saved traces to '" << filename1 << "' and '" << filename2 << "'";
  }
  catch(...)
  {
    mCRL2log(log::info) << ", but its traces could not be saved to '" << filename1 << "' and '" << filename2 << "'";
  }
}

// Facility for constructing a trace to a given state.
template <typename Explorer>
class trace_constructor
{
  protected:
    Explorer& m_explorer;
    std::map<lps::state, lps::state> m_backpointers;

    // Finds a transition s0 --a--> s1, and returns a.
    lps::multi_action find_action(const lps::state& s0, 
                                  const lps::state& s1)
    {
      if constexpr (Explorer::is_stochastic)
      {
        for (const std::pair<lps::multi_action, lps::stochastic_state>& t: m_explorer.generate_transitions(s0))
        {
          for (const lps::state& s: t.second.states)
          {
            if (s == s1)
            {
              return t.first;
            }
          }
        }
      }
      else
      {
        for (const std::pair<lps::multi_action, lps::state>& t: m_explorer.generate_transitions(s0))
        {
          if (t.second == s1)
          {
            return t.first;
          }
        }
      }
      throw mcrl2::runtime_error("no transition found in find_action");
    }

  public:
    explicit trace_constructor(Explorer& explorer_)
      : m_explorer(explorer_)
    {}

    // Constructs a trace ending in s, using the backpointers map.
    class trace construct_trace(const lps::state& s)
    {
      std::deque<lps::state> states{ s };
      std::deque<lps::multi_action> actions;
      while (true)
      {
        const lps::state& s1 = states.front();
        auto i = m_backpointers.find(s1);
        if (i == m_backpointers.end())
        {
          break;
        }
        const lps::state& s0 = i->second;
        states.push_front(s0);
        actions.push_front(find_action(s0, s1));
      }

      class trace tr;
      for (std::size_t i = 0; i < actions.size(); i++)
      {
        tr.set_state(states[i]);
        tr.add_action(actions[i]);
      }
      tr.set_state(states.back());
      return tr;
    }

    // Adds a back pointer for the given edge
    void add_edge(const lps::state& s0, const lps::state& s1)
    {
      m_backpointers[s1] = s0;
    }

    void clear()
    {
      m_backpointers.clear();
    }

    // Providing access to the explorer should perhaps be avoided.
    Explorer& explorer()
    {
      return m_explorer;
    }
};

template <typename Explorer>
class action_detector
{
  protected:
    const std::set<core::identifier_string>& trace_actions;
    const std::set<lps::multi_action>& trace_multiactions;
    trace_constructor<Explorer>& m_trace_constructor;
    const std::string& filename_prefix;
    std::vector<bool> summand_matches;
    std::size_t m_trace_count = 0;
    std::size_t m_max_trace_count;

    bool match_action(const lps::action_summand& summand) const
    {
      using utilities::detail::contains;
      for (const process::action& a: summand.multi_action().actions())
      {
        if (contains(trace_actions, a.label().name()))
        {
          return true;
        }
      }
      return false;
    }

    bool match_summand(std::size_t i) const
    {
      return summand_matches[i];
    }

    std::string create_filename(const lps::multi_action& a)
    {
      using utilities::detail::contains;
      std::string filename = filename_prefix + "_act_" + std::to_string(m_trace_count++);
      if (utilities::detail::contains(trace_multiactions, a))
      {
        filename += "_" + lps::pp(a);
      }
      for (const process::action& a_i: a.actions())
      {
        if (utilities::detail::contains(trace_actions, a_i.label().name()))
        {
          filename += "_" + core::pp(a_i.label().name());
        }
      }
      filename = filename + ".trc";
      return filename;
    }

  public:
    template <typename Specification>
    action_detector(
      const Specification& lpsspec,
      trace_constructor<Explorer>& trace_constructor_,
      const std::set<core::identifier_string>& trace_actions_,
      const std::set<lps::multi_action>& trace_multiactions_,
      const std::string& filename_prefix_,
      std::size_t max_trace_count
    )
      : trace_actions(trace_actions_),
        trace_multiactions(trace_multiactions_),
        m_trace_constructor(trace_constructor_),
        filename_prefix(filename_prefix_),
        m_max_trace_count(max_trace_count)
    {
      using utilities::detail::contains;
      const auto& summands = lpsspec.process().action_summands();
      summand_matches.reserve(summands.size());
      for (const auto& summand: summands)
      {
        summand_matches.push_back(match_action(summand));
      }
    }

    bool detect_action(const lps::state& s0, std::size_t s0_index, const lps::multi_action& a, const lps::state& s1, std::size_t summand_index)
    {
      using utilities::detail::contains;
      if (!match_summand(summand_index))
      {
        return false;
      }
      bool result = false;

      mCRL2log(log::info) << "Action '" + lps::pp(a) + "' found (state index: " + std::to_string(s0_index) + ")";
      if (m_trace_count < m_max_trace_count)
      {
        class trace tr = m_trace_constructor.construct_trace(s0);
        tr.add_action(a);
        tr.set_state(s1);
        std::string filename = create_filename(a);
        save_trace(tr, filename);
        result = true;
      }
      mCRL2log(log::info) << ".\n";
      if (m_max_trace_count > 0 && m_trace_count >= m_max_trace_count)
      {
        m_trace_constructor.explorer().abort();
      }
      return result;
    }
};

template <typename Explorer>
class deadlock_detector
{
  protected:
    trace_constructor<Explorer>& m_trace_constructor;
    const std::string& filename_prefix;
    std::size_t m_trace_count = 0;
    std::size_t m_max_trace_count;

  public:
    deadlock_detector(
      trace_constructor<Explorer>& trace_constructor_,
      const std::string& filename_prefix_,
      std::size_t max_trace_count
    )
      : m_trace_constructor(trace_constructor_),
        filename_prefix(filename_prefix_),
        m_max_trace_count(max_trace_count)
    {}

    void detect_deadlock(const lps::state& s, std::size_t s_index)
    {
      mCRL2log(log::info) << "Deadlock found (state index: " + std::to_string(s_index) + ")";
      if (m_trace_count < m_max_trace_count)
      {
        class trace tr = m_trace_constructor.construct_trace(s);
        std::string filename = filename_prefix + "_dlk_" + std::to_string(m_trace_count++) + ".trc";
        save_trace(tr, filename);
      }
      if (m_max_trace_count > 0 && m_trace_count >= m_max_trace_count)
      {
        m_trace_constructor.explorer().abort();
      }
      mCRL2log(log::info) << ".\n";
    }
};

template <typename Explorer>
class nondeterminism_detector
{
  protected:
    trace_constructor<Explorer>& m_trace_constructor;
    const std::string& filename_prefix;
    std::vector<std::map<lps::multi_action, lps::state> > m_transitions_vec;
    std::size_t m_trace_count = 0;
    std::size_t m_max_trace_count;

  public:
    nondeterminism_detector(
      trace_constructor<Explorer>& trace_constructor_,
      const std::string& filename_prefix_,
      const std::size_t number_of_threads,
      std::size_t max_trace_count = 0
    )
      : m_trace_constructor(trace_constructor_),
        filename_prefix(filename_prefix_),
        m_transitions_vec(number_of_threads+1),  //Threads are number from 1 to n. 
        m_max_trace_count(max_trace_count)
    {
      assert(number_of_threads>0);
    }

    void start_state(std::size_t thread_index)
    {
      assert(thread_index<m_transitions_vec.size());
      m_transitions_vec[thread_index].clear();
    }

    bool detect_nondeterminism(const lps::state& s0, std::size_t s0_index, const lps::multi_action& a, const lps::state& s1, std::size_t thread_index)
    {
      bool result = false;
      assert(thread_index<m_transitions_vec.size());
      auto i = m_transitions_vec[thread_index].find(a);
      if (i == m_transitions_vec[thread_index].end())
      {
        m_transitions_vec[thread_index].insert(std::make_pair(a, s1));
      }
      else if (i->second != s1) // nondeterminism detected
      {
        mCRL2log(log::info) << "Nondeterministic state found (state index: " + std::to_string(s0_index) + ")";
        if (m_trace_count < m_max_trace_count)
        {
          class trace tr = m_trace_constructor.construct_trace(s0);
          tr.add_action(a);
          tr.set_state(s1);
          std::string filename = filename_prefix + "_nondeterministic_" + std::to_string(m_trace_count++) + ".trc";
          save_trace(tr, filename);
          result = true;
        }
        mCRL2log(log::info) << ".\n";
        if (m_max_trace_count > 0 && m_trace_count >= m_max_trace_count)
        {
          m_trace_constructor.explorer().abort();
        }
      }
      return result;
    }
};

template <typename Explorer>
class divergence_detector
{
  public:
    using state_type = typename Explorer::state_type;
    using state_index_type = typename Explorer::state_index_type;

    // data type for storing the last discovered states
    using last_discovered_type = std::conditional_t<Explorer::is_stochastic, std::forward_list<lps::state>, lps::state>;

  protected:
    Explorer& explorer;
    const std::string& filename_prefix;
    trace_constructor<Explorer> m_local_trace_constructor;
    utilities::unordered_map<lps::state, std::size_t> m_divergent_states;
    std::vector<lps::explorer_summand> m_regular_summands;
    std::vector<lps::explorer_summand> m_confluent_summands;
    std::size_t m_trace_count = 0;
    std::size_t m_max_trace_count;
    std::mutex divergence_detector_mutex;  // As it stands the divergence detector is sequential. 

  public:
    divergence_detector(
      Explorer& explorer_,
      const std::set<core::identifier_string>& actions,
      const std::string& filename_prefix_,
      std::size_t max_trace_count
    )
      : explorer(explorer_),
        filename_prefix(filename_prefix_),
        m_local_trace_constructor(explorer),
        m_max_trace_count(max_trace_count)
    {
      using utilities::detail::contains;

      auto is_hidden = [&](const lps::explorer_summand& summand)
      {
        for (const process::action& a: summand.multi_action.actions())
        {
          if (!contains(actions, a.label().name()))
          {
            return false;
          }
        }
        return true;
      };

      for (const lps::explorer_summand& summand: explorer.regular_summands())
      {
        if (is_hidden(summand))
        {
          m_regular_summands.push_back(summand);
        }
      }

      for (const lps::explorer_summand& summand: explorer.confluent_summands())
      {
        if (is_hidden(summand))
        {
          m_confluent_summands.push_back(summand);
        }
      }
    }

    // Returns true if a trace was saved.
    bool detect_divergence(const lps::state& s, std::size_t s_index, trace_constructor<Explorer>& global_trace_constructor, bool dfs_recursive = false)
    {
      using utilities::detail::contains;

      bool result = false;
      divergence_detector_mutex.lock();
      m_local_trace_constructor.clear();

      auto q = m_divergent_states.find(s);
      if (q != m_divergent_states.end())
      {
        std::string message = "Divergent state found (state index: " + std::to_string(s_index) +
                              "), reachable from divergent state with index " + std::to_string(q->second);
        mCRL2log(log::info) << message << ".\n";
        m_divergent_states.erase(q);
        divergence_detector_mutex.unlock();
        return false;
      }

      std::unordered_set<lps::state> discovered;
      data::data_expression_list process_parameter_undo = explorer.process_parameter_values();

      if (dfs_recursive)
      {
        std::unordered_set<lps::state> gray;
        explorer.generate_state_space_dfs_recursive(
          s,
          gray,
          discovered,
          m_regular_summands,
          m_confluent_summands,
          utilities::skip(), // discover_state
          utilities::skip(), // examine_transition
          utilities::skip(), // tree_edge

          // back_edge
          [&](const lps::state& s0, const lps::multi_action& a, const state_type& s1) {
            mCRL2log(log::info) << "Divergent state found (state index: " + std::to_string(s_index) + ")";
            if (m_trace_count < m_max_trace_count)
            {
              class trace tr = global_trace_constructor.construct_trace(s);
              class trace tr_loop = m_local_trace_constructor.construct_trace(s0);
              for (const lps::state& u: tr_loop.states())
              {
                m_divergent_states[u] = s_index;
              }
              tr_loop.add_action(a);
              tr_loop.set_state(first_state(s1));
              std::string filename = filename_prefix + "_divergence_" + std::to_string(m_trace_count) + ".trc";
              std::string loop_filename = filename_prefix + "_divergence_loop" + std::to_string(m_trace_count++) + ".trc";
              save_traces(tr, filename, tr_loop, loop_filename);
              result = true;
            }
            mCRL2log(log::info) << ".\n";
            //--- Workaround for Visual Studio 2019 ---//
            // explorer.abort();
            static_cast<lps::abortable&>(explorer).abort();
          }
        );
      }
      else
      {
        explorer.generate_state_space_dfs_iterative(
          s,
          discovered,
          m_regular_summands,
          m_confluent_summands,
          utilities::skip(), // discover_state
          utilities::skip(), // examine_transition
          utilities::skip(), // tree_edge

          // back_edge
          [&](const lps::state& s0, const lps::multi_action& a, const state_type& s1) {
            mCRL2log(log::info) << "Divergent state found (state index: " + std::to_string(s_index) + ")";
            if (m_trace_count < m_max_trace_count)
            {
              class trace tr = global_trace_constructor.construct_trace(s);
              class trace tr_loop = m_local_trace_constructor.construct_trace(s0);
              for (const lps::state& u: tr_loop.states())
              {
                m_divergent_states[u] = s_index;
              }
              tr_loop.add_action(a);
              tr_loop.set_state(first_state(s1));
              std::string filename = filename_prefix + "_divergence_" + std::to_string(m_trace_count) + ".trc";
              std::string loop_filename = filename_prefix + "_divergence_loop" + std::to_string(m_trace_count++) + ".trc";
              save_traces(tr, filename, tr_loop, loop_filename);
              result = true;
            }
            mCRL2log(log::info) << ".\n";
            //--- Workaround for Visual Studio 2019 ---//
            // explorer.abort();
            static_cast<lps::abortable&>(explorer).abort();
          }
        );
      }
      explorer.set_process_parameter_values(process_parameter_undo);
      if (m_max_trace_count > 0 && m_trace_count >= m_max_trace_count)
      {
        explorer.abort();
      }
      divergence_detector_mutex.unlock();
      return result;
    }
};

class progress_monitor
{
  protected:
    std::size_t level = 1;    // the current exploration level
    std::size_t level_up = 1; // when count reaches level_up, the level is increased
    std::atomic<std::size_t> count = 0;
    std::atomic<std::size_t> transition_count = 0;

    std::size_t last_state_count = 0;
    std::size_t last_transition_count = 0;
    std::atomic<time_t> last_log_time = time(nullptr) - 1;

    lps::exploration_strategy search_strategy;

  public:
    explicit progress_monitor(lps::exploration_strategy search_strategy_)
      : search_strategy(search_strategy_)
    {}

    void examine_transition()
    {
      transition_count++;
    }

    void finish_state(std::size_t state_count, std::size_t todo_list_size, std::size_t number_of_threads)
    {
      time_t new_log_time = 0;

      static std::mutex exclusive_print_mutex;
      if (search_strategy == lps::es_breadth)
      {
        ++count;
        if (number_of_threads == 1 && count == level_up) 
        {
          exclusive_print_mutex.lock();
          mCRL2log(log::debug) << "Number of states at level " << level << " is " << state_count - last_state_count << "\n";
          level++;
          level_up = count + todo_list_size;
          last_state_count = state_count;
          last_transition_count = transition_count;
          exclusive_print_mutex.unlock();
        }

        if (time(&new_log_time) > last_log_time.load(std::memory_order_relaxed))
        {
          exclusive_print_mutex.lock();

          last_log_time = new_log_time;
          std::size_t lvl_states = state_count - last_state_count;
          std::size_t lvl_transitions = transition_count - last_transition_count;
          if (number_of_threads>1) // Levels have no meaning with multiple threads. 
          {
            mCRL2log(log::status) << std::fixed << std::setprecision(2)
                                  << state_count << "st, " << transition_count << "tr"
                                  << ", explored " << 100.0 * ((float) count / state_count)
                                  << "%.\n";
          }
          else
          {
            mCRL2log(log::status) << std::fixed << std::setprecision(2)
                                  << state_count << "st, " << transition_count << "tr"
                                  << ", explored " << 100.0 * ((float) count / state_count)
                                  << "%. Last level: " << level << ", " << lvl_states << "st, " 
                                  << lvl_transitions << "tr.\n";
          }
          exclusive_print_mutex.unlock();
        }
      }
      else
      {
        count++;
        if (time(&new_log_time) > last_log_time.load(std::memory_order_relaxed))
        {
          exclusive_print_mutex.lock();
          last_log_time = new_log_time;
          mCRL2log(log::status) << "monitor: currently explored "
                            << count << " state" << ((count==1)?"":"s")
                            << " and " << transition_count << " transition" << ((transition_count==1)?".":"s.")
                            << std::endl;
          exclusive_print_mutex.unlock();
        }
      }
    }

    void finish_exploration(std::size_t state_count, std::size_t number_of_threads)
    {
      if (search_strategy == lps::es_breadth)
      {
        mCRL2log(log::verbose) << "Done with state space generation (";
        if (number_of_threads==1)
        {
          mCRL2log(log::verbose) << level-1 << " level" << ((level==2)?"":"s") << ", ";
        }
        mCRL2log(log::verbose) << state_count << " state" << ((state_count == 1)?"":"s")
                               << " and " << transition_count << " transition" << ((transition_count==1)?"":"s") << ")" << std::endl;
      }
      else
      {
        mCRL2log(log::verbose) << "Done with state space generation ("
                          << state_count << " state" << ((state_count == 1)?"":"s")
                          << " and " << transition_count << " transition" << ((transition_count==1)?"":"s") << ")" << std::endl;
      }
    }
};

} // namespace detail

template <bool Stochastic, bool Timed, typename Specification>
struct state_space_generator
{
  using explorer_type = lps::explorer<Stochastic, Timed, Specification>;
  using state_type = typename explorer_type::state_type;

  const lps::explorer_options& options;
  lps::explorer<Stochastic, Timed, Specification> explorer;
  detail::trace_constructor<explorer_type> m_trace_constructor;

  detail::action_detector<explorer_type> m_action_detector;
  detail::deadlock_detector<explorer_type> m_deadlock_detector;
  detail::nondeterminism_detector<explorer_type> m_nondeterminism_detector;
  std::unique_ptr<detail::divergence_detector<explorer_type>> m_divergence_detector;
  detail::progress_monitor m_progress_monitor;

  state_space_generator(const Specification& lpsspec, const lps::explorer_options& options_)
    : options(options_),
      explorer(lpsspec, options_),
      m_trace_constructor(explorer),
      m_action_detector(lpsspec, m_trace_constructor, options.trace_actions, options.trace_multiactions, options.trace_prefix, options.max_traces),
      m_deadlock_detector(m_trace_constructor, options.trace_prefix, options.max_traces),
      m_nondeterminism_detector(m_trace_constructor, options.trace_prefix, options.number_of_threads, options.max_traces),
      m_progress_monitor(options.search_strategy)
  {
    if (options.detect_divergence)
    {
      m_divergence_detector = 
               std::unique_ptr<detail::divergence_detector<explorer_type>>(
                        new detail::divergence_detector<explorer_type>(explorer, 
                                                                       options.actions_internal_for_divergencies, 
                                                                       options.trace_prefix, 
                                                                       options.max_traces));
    }
  }

  bool max_states_exceeded(const std::size_t thread_index)
  {
    return explorer.state_map().size(thread_index) >= options.max_states;
  }

  struct aligned_bool 
  {
    alignas(64) size_t m_bool;
  };

  // Explore the specification passed via the constructor, and put the results in builder.
  template <typename LTSBuilder>
  bool explore(LTSBuilder& builder)
  {
    std::vector<aligned_bool> has_outgoing_transitions(options.number_of_threads+1); // thread indices start at 1. 
    const lps::state* source = nullptr;

    try
    {
      explorer.generate_state_space(
        false,

        // discover_state
        [&](const std::size_t thread_index, const lps::state& s, std::size_t s_index)
        {
          if (options.generate_traces && source)
          {
            m_trace_constructor.add_edge(*source, s);
          }
          if (options.detect_divergence)
          {
            // TODO: support divergence checks for stochastic specifications
            if constexpr (!Stochastic)
            {
              m_divergence_detector->detect_divergence(s, s_index, m_trace_constructor, options.dfs_recursive);
            }
          }
          // if (explorer.state_map().size() >= options.max_states)
          //--- Workaround for Visual Studio 2019 ---//
          if (max_states_exceeded(thread_index))
          {
            static bool not_reported_yet=true;
            if (not_reported_yet)
            {
              not_reported_yet=false;
              mCRL2log(log::verbose) << "Explored the maximum number (" << options.max_states << ") of states, terminating." << std::endl;
            }
            //--- Workaround for Visual Studio 2019 ---//
            // explorer.abort();
            static_cast<lps::abortable&>(explorer).abort();
          }
        },

        // examine_transition
        [&](const std::size_t thread_index, const std::size_t number_of_threads, 
            const lps::state& s0, std::size_t s0_index, const lps::multi_action& a, 
            const auto& s1, const auto& s1_index, std::size_t summand_index)
        {
          if constexpr (Stochastic)
          {
            builder.add_transition(s0_index, a, s1_index, s1.probabilities, number_of_threads);
          }
          else
          {
            builder.add_transition(s0_index, a, s1_index, number_of_threads);
          }
          assert(thread_index<has_outgoing_transitions.size());
          has_outgoing_transitions[thread_index].m_bool = true;
          if (options.detect_action)
          {
            m_action_detector.detect_action(s0, s0_index, a, first_state(s1), summand_index);
          }
          if (options.detect_nondeterminism)
          {
            m_nondeterminism_detector.detect_nondeterminism(s0, s0_index, a, first_state(s1), thread_index);
          }
          if (!options.suppress_progress_messages)
          {
            m_progress_monitor.examine_transition();
          }
        },

        // start_state
        [&](const std::size_t thread_index, const lps::state& s, std::size_t /* s_index */)
        {
          if (options.number_of_threads == 1) {
            source = &s;
          }

          assert(thread_index<has_outgoing_transitions.size());
          has_outgoing_transitions[thread_index].m_bool = false;
          if (options.detect_nondeterminism)
          {
            m_nondeterminism_detector.start_state(thread_index);
          }
        },

        // finish_state
        [&](const std::size_t thread_index, const std::size_t number_of_threads, 
            const lps::state& s, std::size_t s_index, std::size_t todo_list_size)
        {
          assert(thread_index<has_outgoing_transitions.size());
          if (options.detect_deadlock && !has_outgoing_transitions[thread_index].m_bool)
          {
            m_deadlock_detector.detect_deadlock(s, s_index);
          }
          if (!options.suppress_progress_messages)
          {
            m_progress_monitor.finish_state(explorer.state_map().size(thread_index), todo_list_size, number_of_threads);
          }
        },

        // discover_initial_state
        [&](const lps::stochastic_state& s, const std::list<std::size_t>& s_index)
        {
          if constexpr (Stochastic)
          {
            builder.set_initial_state(s_index, s.probabilities);
          }
        }
      );
      m_progress_monitor.finish_exploration(explorer.state_map().size(), options.number_of_threads);
      builder.finalize(explorer.state_map(), Timed);
    }
    catch (const data::enumerator_error& e)
    {
      mCRL2log(log::error) << "Error while exploring state space: " << e.what() << ".\n";
      if (options.save_error_trace)
      {
        const lps::state& s = *source;
        class trace tr = m_trace_constructor.construct_trace(s);
        std::string filename = options.trace_prefix + "_error.trc";
        detail::save_trace(tr, filename);
      }
      return false;
    }

    return true;
  }
};

} // namespace mcrl2::lts

#endif // MCRL2_LTS_STATE_SPACE_GENERATOR_H
