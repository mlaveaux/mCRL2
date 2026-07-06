// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/lts/stochastic_lts_builder.h
/// \brief add your file description here.

#ifndef MCRL2_LTS_STOCHASTIC_LTS_BUILDER_H
#define MCRL2_LTS_STOCHASTIC_LTS_BUILDER_H

#include "mcrl2/lts/lts_builder.h"

namespace mcrl2::lts {

struct stochastic_lts_builder
{
  using indexed_set_for_states_type = atermpp::indexed_set<lps::state, mcrl2::utilities::detail::GlobalThreadSafe>;
  // All LTS classes use integers to represent actions in transitions. A mapping from actions to integers
  // is needed to avoid duplicates.
  utilities::unordered_map_large<lps::multi_action, std::size_t> m_actions;

  stochastic_lts_builder()
  {
    lps::multi_action tau(process::action_list(), data::undefined_real());
    m_actions.emplace(std::make_pair(tau, m_actions.size()));
  }

  std::size_t add_action(const lps::multi_action& a)
  {
    auto i = m_actions.find(a);
    if (i == m_actions.end())
    {
      i = m_actions.emplace(std::make_pair(a, m_actions.size())).first;
    }
    return i->second;
  }

  // Set the initial (stochastic) state of the LTS
  virtual void set_initial_state(const std::list<std::size_t>& targets, const std::vector<data::data_expression>& probabilities) = 0;

  // Add a transition to the LTS
  virtual void add_transition(std::size_t from,
    const lps::multi_action& a,
    const std::list<std::size_t>& targets,
    const std::vector<data::data_expression>& probabilities,
    std::size_t number_of_threads = 1)
    = 0;

  // Add actions and states to the LTS
  virtual void finalize(const indexed_set_for_states_type& state_map, bool timed) = 0;

  // Save the LTS to a file
  virtual void save(const std::string& filename) = 0;

  virtual ~stochastic_lts_builder() = default;
};

class stochastic_lts_none_builder: public stochastic_lts_builder
{
  public:
    void set_initial_state(const std::list<std::size_t>& /* to */, const std::vector<data::data_expression>& /* probabilities */) override
    {}

    void add_transition(std::size_t /* from */, const lps::multi_action& /* a */, const std::list<std::size_t>& /* targets */, const std::vector<data::data_expression>& /* probabilities */, const std::size_t /* number_of_threads */) override
    {}

    void finalize(const indexed_set_for_states_type& /* state_map */, bool /* timed */) override
    {}

    void save(const std::string& /* filename */) override
    {}
};

class stochastic_lts_aut_builder: public stochastic_lts_builder
{
  protected:
    struct stochastic_state
    {
      std::list<std::size_t> targets;
      std::vector<data::data_expression> probabilities;

      stochastic_state() = default;

      stochastic_state(std::list<std::size_t>  targets_, std::vector<data::data_expression>  probabilities_)
        : targets(std::move(targets_)), probabilities(std::move(probabilities_))
      {}

      void save_to_aut(std::ostream& out) const
      {
        auto j = targets.begin();
        out << *j++;
        if (targets.size() > 1)
        {
          for (auto i = probabilities.begin(); j != targets.end(); ++i, ++j)
          {
            out << " " << lps::print_probability(*i) << " " << *j;
          }
        }
      }
    };

    struct transition
    {
      std::size_t from;
      std::size_t label;
      std::size_t to;

      transition(std::size_t from_, std::size_t label_, std::size_t to_)
        : from(from_), label(label_), to(to_)
      {}

      bool operator<(const transition& other) const
      {
        return std::tie(from, label, to) < std::tie(other.from, other.label, other.to);
      }
    };

    std::vector<stochastic_state> m_stochastic_states;
    std::vector<transition> m_transitions;
    std::size_t m_number_of_states = 0;
    std::mutex m_exclusive_transition_access;

  public:
    stochastic_lts_aut_builder() = default;

    // Set the initial (stochastic) state of the LTS
    void set_initial_state(const std::list<std::size_t>& targets, const std::vector<data::data_expression>& probabilities) override
    {
      m_stochastic_states.emplace_back(targets, probabilities);
    }

    // Add a transition to the LTS
    void add_transition(std::size_t from, const lps::multi_action& a, const std::list<std::size_t>& targets, const std::vector<data::data_expression>& probabilities, const std::size_t number_of_threads) override
    {
      std::unique_lock<std::mutex> guard(m_exclusive_transition_access, std::defer_lock);
      if (mcrl2::utilities::detail::GlobalThreadSafe && number_of_threads > 1)
      {
        guard.lock();
      }
      std::size_t to = m_stochastic_states.size();
      std::size_t label = add_action(a);
      m_stochastic_states.emplace_back(targets, probabilities);
      m_transitions.emplace_back(from, label, to);
    }

    // Add actions and states to the LTS
    void finalize(const indexed_set_for_states_type& state_map, bool /* timed */) override
    {
      m_number_of_states = state_map.size();
    }

    void save(std::ostream& out) const
    {
      std::vector<lps::multi_action> actions{ m_actions.size() };
      for (const auto& p: m_actions)
      {
        actions[p.second] = p.first;
      }

      out << "des (";
      m_stochastic_states[0].save_to_aut(out);
      out << "," << m_transitions.size() << "," << m_number_of_states << ")" << "\n";

      for (const transition& t: m_transitions)
      {
        out << "(" << t.from << ",\"" << lps::pp(actions[t.label]) << "\",";
        m_stochastic_states[t.to].save_to_aut(out);
        out << ")" << "\n";
      }
    }

    // Save the LTS to a file
    void save(const std::string& filename) override
    {
      if (filename.empty())
      {
        save(std::cout);
      }
      else
      {
        std::ofstream out(filename.c_str());
        if (!out.is_open())
        {
          throw mcrl2::runtime_error("cannot create .aut file '" + filename + ".");
        }
        save(out);
        out.close();
      }
    }
};

class stochastic_lts_lts_builder: public stochastic_lts_builder
{
  protected:
    probabilistic_lts_lts_t m_lts;
    bool m_discard_state_labels = false;
    probabilistic_state<std::size_t, lps::probabilistic_data_expression> m_initial_state;
    std::mutex m_exclusive_transition_access;

  public:
    stochastic_lts_lts_builder(
      const data::data_specification& dataspec,
      const process::action_label_list& action_labels,
      const data::variable_list& process_parameters,
      bool discard_state_labels = false
    )
      : m_discard_state_labels(discard_state_labels)
    {
      m_lts.set_data(dataspec);
      m_lts.set_process_parameters(process_parameters);
      m_lts.set_action_label_declarations(action_labels);
    }

    static probabilistic_state<std::size_t, lps::probabilistic_data_expression> 
             make_probabilistic_state(const std::list<std::size_t>& targets, const std::vector<data::data_expression>& probabilities)
    {
      probabilistic_state<std::size_t, lps::probabilistic_data_expression> result;
      assert(targets.size()>0);
      if (targets.size()==1)
      {
        result.set(targets.front());
      }
      else 
      {
        std::list<std::size_t>::const_iterator ti = targets.begin();
        std::vector<data::data_expression>::const_iterator pi = probabilities.begin();
        for (; ti != targets.end(); ++pi, ++ti)
        {
          result.add(*ti, lps::probabilistic_data_expression(*pi));
        }
      }
      return result;
    }

    // Set the initial (stochastic) state of the LTS
    void set_initial_state(const std::list<std::size_t>& targets, const std::vector<data::data_expression>& probabilities) override
    {
      m_initial_state = make_probabilistic_state(targets, probabilities);
    }

    // Add a transition to the LTS
    void add_transition(std::size_t from, const lps::multi_action& a, const std::list<std::size_t>& targets, const std::vector<data::data_expression>& probabilities, const std::size_t number_of_threads) override
    {
      std::unique_lock<std::mutex> guard(m_exclusive_transition_access, std::defer_lock);
      if (mcrl2::utilities::detail::GlobalThreadSafe && number_of_threads > 1)
      {
        guard.lock();
      }
      auto s1 = make_probabilistic_state(targets, probabilities);
      std::size_t label = add_action(a);
      std::size_t to = m_lts.add_probabilistic_state(s1);
      m_lts.add_transition(transition(from, label, to));
    }

    // Add actions and states to the LTS
    void finalize(const indexed_set_for_states_type& state_map, bool timed) override
    {
      // add actions
      m_lts.set_num_action_labels(m_actions.size());
      for (const auto& p: m_actions)
      {
        m_lts.set_action_label(p.second, action_label_lts(lps::multi_action(p.first.actions(), p.first.time())));
      }

      // add state labels
      if (!m_discard_state_labels)
      {
        std::size_t n = state_map.size();
        std::vector<state_label_lts> state_labels(n);
        for (std::size_t i = 0; i < n; i++)
        {
          if (timed)
          {
            state_labels[i] = state_label_lts(remove_time_stamp(state_map[i]));
          }
          else
          {
            state_labels[i] = state_label_lts(state_map[i]);
          }
        }
        m_lts.state_labels() = std::move(state_labels);
      }

      m_lts.set_num_states(state_map.size(), true);
      m_lts.set_initial_probabilistic_state(m_initial_state); // This can't be done at the start :-(
    }

    // Save the LTS to a file
    void save(const std::string& filename) override
    {
      m_lts.save(filename);
    }
};

class stochastic_lts_fsm_builder: public stochastic_lts_lts_builder
{
  public:
    using super = stochastic_lts_lts_builder;
    stochastic_lts_fsm_builder(const data::data_specification& dataspec, const process::action_label_list& action_labels, const data::variable_list& process_parameters)
      : super(dataspec, action_labels, process_parameters)
    { }

    void save(const std::string& filename) override
    {
      probabilistic_lts_fsm_t fsm;
      detail::lts_convert(m_lts, fsm);
      fsm.save(filename);
    }
};

inline
std::unique_ptr<stochastic_lts_builder> create_stochastic_lts_builder(const lps::stochastic_specification& lpsspec, const lps::explorer_options& options, lts_type output_format)
{
  switch (output_format)
  {
    case lts_aut: return std::make_unique<stochastic_lts_aut_builder>();
    case lts_lts: return std::make_unique<stochastic_lts_lts_builder>(lpsspec.data(), lpsspec.action_labels(), lpsspec.process().process_parameters(), options.discard_lts_state_labels);
    case lts_fsm: return std::make_unique<stochastic_lts_fsm_builder>(lpsspec.data(), lpsspec.action_labels(), lpsspec.process().process_parameters());
    default: return std::make_unique<stochastic_lts_none_builder>();
  }
}

} // namespace mcrl2::lts



#endif // MCRL2_LTS_STOCHASTIC_LTS_BUILDER_H
