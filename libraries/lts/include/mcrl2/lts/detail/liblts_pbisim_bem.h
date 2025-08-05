// Author(s): Hector Joao Rivera Verduzco
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file lts/detail/liblts_pbisim_bem.h

#ifndef MCRL2_LTS_DETAIL_LIBLTS_PBISIM_BEM_H
#define MCRL2_LTS_DETAIL_LIBLTS_PBISIM_BEM_H

#include "mcrl2/lts/detail/liblts_plts_merge.h"

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/execution_timer.cxx"
#else  
  import utilities;
#endif

namespace mcrl2::lts::detail
{

template < class LTS_TYPE>
class prob_bisim_partitioner_bem
{

  public:
    /** \brief Creates a probabilistic bisimulation partitioner for an PLTS.
     *  \details This bisimulation partitioner applies the algorithm defined in C. Baier, B. Engelen and M. Majster-Cederbaum. 
     *   Deciding Bisimilarity and Similarity for Probabilistic Processes. In Journal of Computer and System Sciences 60, 187-237 (2000)
     */
    prob_bisim_partitioner_bem(LTS_TYPE& l, utilities::execution_timer& timer):aut(l)
    {
      mCRL2log(log::verbose) << "Probabilistic bisimulation partitioner created for "
                  << l.num_states() << " states and " <<
                  l.num_transitions() << " transitions\n";
      timer.start("bisimulation_reduce (bem)");
      create_initial_partition();
      refine_partition_until_it_becomes_stable();
      timer.finish("bisimulation_reduce (bem)");
      postprocessing_stage();
    }

    /** \brief Gives the number of bisimulation equivalence classes of the LTS.
    *  \return The number of bisimulation equivalence classes of the LTS.
    */
    std::size_t num_eq_classes() const
    {
      return state_partition.size();
    }

    /** \brief Gives the bisimulation equivalence class number of a state.
     *  \param[in] s A state number.
     *  \return The number of the bisimulation equivalence class to which \e s belongs. */
    std::size_t get_eq_class(const std::size_t s) const
    {
      assert(s<block_index_of_a_state.size());
      return block_index_of_a_state[s]; // The block index is the state number of the block.
    }

    /** \brief Gives the bisimulation equivalence step class number of a probabilistic state.
    *  \param[in] s A probabilistic state number.
    *  \return The number of the step class to which s belongs. */
    std::size_t get_eq_step_class(const std::size_t d) const
    {
      assert(d < step_class_index_of_a_distribution.size());
      assert(step_class_index_of_a_distribution[d] < step_classes.size());
      assert(step_classes[step_class_index_of_a_distribution[d]].equivalent_step_class < step_classes.size());
      return step_classes[step_class_index_of_a_distribution[d]].equivalent_step_class;
    }

    /** \brief Replaces the transition relation of the current lts by the transitions
    *         of the bisimulation reduced transition system.
    * \pre The bisimulation equivalence classes have been computed. */
    void replace_transitions()
    {
      std::set < transition > resulting_transitions;

      const std::vector<transition>& trans = aut.get_transitions();
      for (const transition& t : trans)
      {
        resulting_transitions.insert(
          transition(
            get_eq_class(t.from()),
            t.label(),
            get_eq_step_class(t.to())));
      }
      // Remove the old transitions
      aut.clear_transitions();

      // Copy the transitions from the set into the transition system.
      for (const transition& t : resulting_transitions)
      {
        aut.add_transition(t);
      }
    }

    /** \brief Replaces the probabilistic states of the current lts by the probabilistic
    *         states of the bisimulation reduced transition system.
    * \pre The bisimulation step classes have been computed. */
    void replace_probabilistic_states()
    {
      std::vector<typename LTS_TYPE::probabilistic_state_t> new_probabilistic_states;

      // get the equivalent probabilistic state of each probabilistic block and add it to aut
      for (step_class_type& sc : equivalent_step_classes)
      {
        typename LTS_TYPE::probabilistic_state_t equivalent_ps = calculate_equivalent_probabilistic_state(sc);
        new_probabilistic_states.push_back(equivalent_ps);
      }

      /* Remove old probabilistic states */
      aut.clear_probabilistic_states();

      // Add new prob states to aut
      for (const typename LTS_TYPE::probabilistic_state_t& new_ps : new_probabilistic_states)
      {
        aut.add_probabilistic_state(new_ps);
      }

      typename LTS_TYPE::probabilistic_state_t old_initial_prob_state = aut.initial_probabilistic_state();
      typename LTS_TYPE::probabilistic_state_t new_initial_prob_state = calculate_new_probabilistic_state(old_initial_prob_state);
      aut.set_initial_probabilistic_state(new_initial_prob_state);
    }

    /** \brief Returns whether two states are in the same probabilistic bisimulation equivalence class.
    *  \param[in] s A state number.
    *  \param[in] t A state number.
    *  \retval true if \e s and \e t are in the same bisimulation equivalence class;
    *  \retval false otherwise. */
    bool in_same_probabilistic_class(const std::size_t s, const std::size_t t) const
    {
      return get_eq_step_class(s) == get_eq_step_class(t);
    }

  private:

    using block_key_type = std::size_t;
    using step_class_key_type = std::size_t;
    using state_type = std::size_t;
    using label_type = std::size_t;
    using distribution_key_type = std::size_t;
    using probability_fraction_type = typename LTS_TYPE::probabilistic_state_t::probability_t;

    struct distribution_type
    {
      distribution_key_type key = 0UL;
      std::vector<std::list<transition*>> incoming_transitions_per_label; // Incoming transitions organized per label
  };
  
  struct block_type
  {
    block_key_type key = 0UL;
    std::list<state_type> states;        // The states in the block
    bool is_in_new_blocks = false;
  };

  struct step_class_type {
    step_class_key_type key = 0UL;
    label_type action = 0UL;                          // action label of the pair <a,M>.
    std::list<distribution_type*> distributions;      // The distributions in the step class <a,M>.
    std::vector< bool > prev_states;                  // Previous states that can reach step_class <a,M>
    bool is_in_new_step_classes = false;
    std::size_t equivalent_step_class = 0UL;
  };

  std::list<block_type*> state_partition;
  std::list<step_class_type*> step_partition;
  std::vector<distribution_type> distributions;
  std::deque<step_class_type> step_classes;
  std::deque<step_class_type> equivalent_step_classes;
  std::deque<block_type> blocks;
  std::vector<block_key_type> block_index_of_a_state;
  std::vector<step_class_key_type> step_class_index_of_a_distribution;

  LTS_TYPE& aut;

  struct tree_type
  {
    std::list<state_type> states;
    std::size_t count = 0UL;
    tree_type* left = nullptr;
    tree_type* right = nullptr;

    void init_node()
    {
      count = 0;
      left = nullptr;
      right = nullptr;
    }
  };

  /** \brief Creates the initial partition of step classes and blocks.
   *  \detail The blocks are initially partitioned based on the actions that can perform.
   *          The step classes are partitioned based on the action that leads to the probabilistic state */
  void create_initial_partition()
  {
    std::vector< std::vector< std::list<distribution_type*> > > steps; // Representation of transition in 2-d array
    std::vector<transition>& transitions = aut.get_transitions();
    std::vector< std::vector<bool> > distribution_per_step_class;
    
    distribution_per_step_class.resize(aut.num_action_labels());
    for (std::vector<bool>& i : distribution_per_step_class)
    {
      i.resize(aut.num_probabilistic_states());
    }

    //---- Preprocessing stage to transform the PLTS to the data structures suggested by Baier ---- //

    // Construct vector of distributions
    distributions.resize(aut.num_probabilistic_states());
    std::size_t key = 0;
    for (distribution_type& distribution : distributions)
    {
      distribution.key = key;
      distribution.incoming_transitions_per_label.resize(aut.num_action_labels());
      key++;
    }

    // Initialize the Steps bi-dimientional array
    steps.resize(aut.num_states());
    for (std::size_t i = 0; i < steps.size(); i++)
    {
      steps[i].resize(aut.num_action_labels());
    }

    for (transition& t : transitions)
    {
      steps[t.from()][t.label()].push_back(&distributions[t.to()]);
      distributions[t.to()].incoming_transitions_per_label[t.label()].push_back(&t);
    }

    //---- Start the initialization process (page 207. Fig. 10. Baier) ---- //

    // Initially there are as many step classes as lables
    step_classes.resize(aut.num_action_labels());

    for (step_class_type& sc : step_classes)
    {
      sc.prev_states.resize(aut.num_states());
    }

    // create tree structure to group action states based on the outgoing transitions
    std::size_t max_block_size = 0;
    tree_type* max_block = nullptr;
    std::deque<tree_type> tree_nodes;
    std::vector<tree_type*> leaves;
    tree_type v0;
    v0.init_node();
    tree_nodes.push_back(v0);

    // For all s in S do
    for (state_type s = 0; s < aut.num_states(); s++)
    {
      // (1) Create pointer to root of the tree
      tree_type* v = &tree_nodes[0];

      // (2) Construct tree
      for (std::size_t i = 0; i < aut.num_action_labels(); i++)
      {
        step_classes[i].key = i;
        step_classes[i].action = i;
        if (steps[s][i].size() > 0)
        {
          for (distribution_type* d : steps[s][i])
          {
            if (distribution_per_step_class[i][d->key] == false)
            {
              step_classes[i].distributions.push_back(d);
              distribution_per_step_class[i][d->key] = true;
            }
          }
          
          step_classes[i].prev_states[s] = true;

          if (v->left == nullptr)
          {
            // create left node
            tree_type w;
            w.init_node();
            tree_nodes.push_back(w);
            v->left = &tree_nodes.back();

            // add new node to the leaf nodes if it is a leaf
            if (i == aut.num_action_labels() - 1)
            {
              leaves.push_back(v->left);
            }

          }
          v = v->left;
        }
        else
        {
          if (v->right == nullptr)
          {
            // create left node
            tree_type w;
            w.init_node();
            tree_nodes.push_back(w);
            v->right = &tree_nodes.back();

            // add new node to the leaf nodes if it is a leaf
            if (i == aut.num_action_labels() - 1)
            {
              leaves.push_back(v->right);
            }
          }
          v = v->right;
        }
      }

      // (3) Add the state to the leaf and increment its state counter
      v->states.push_back(s);
      v->count++;

      // (4) Keep track of the leave containing more states
      if (v->count > max_block_size)
      {
        max_block = v;
        max_block_size = v->count;
      }
    }

    // insert all states of the leaves to blocks
    blocks.resize(leaves.size());

    key = 0;
    std::size_t larger_key = 0;
    for (tree_type* leaf_ptr : leaves)
    {
      block_type& block = blocks[key];
      block.key = key;
      block.states.splice(block.states.end(), leaf_ptr->states);

      if (leaf_ptr == max_block)
      {
        larger_key = key;
      }

      key++;
    }

    // Add all blocks to the state partition.
    // Initially only the larger block has to be at the end of the list
    for (block_type& b : blocks)
    {
      if (b.key != larger_key)
      {
        // Push the non-larger blocks to the front of the list. This are the so called new blocks
        b.is_in_new_blocks = true;
        state_partition.push_front(&b);
      }
      else
      {
        // push the larger block to the end of the list
        b.is_in_new_blocks = false;
        state_partition.push_back(&b);
      }
    }
    
    // Add all step classes to the step partition.
    // Initially all are new step classes
    for (step_class_type& sc : step_classes)
    {
      if (sc.distributions.size() > 0)
      { 
        step_partition.push_front(&sc);
        sc.is_in_new_step_classes = true;
      }
    }

    block_index_of_a_state.resize(aut.num_states());
    for (const block_type& b : blocks)
    {
      for (const state_type s : b.states)
      {
        block_index_of_a_state[s] = b.key;
      }
    }

  }

  /** \brief Calculates the probability to reach block b from distribution d.
  *   \param[in] d is a probabilistic state (distribution).
  *              b is a block of states.
  *   \return The probability to reach block b. */
  probability_fraction_type probability_to_block(distribution_type& d, block_type& b)
  {
    probability_fraction_type prob_to_block;
    const typename LTS_TYPE::probabilistic_state_t& prob_state = aut.probabilistic_state(d.key);

    /* Check whether the state is in the distribution d and add up the probability*/
    if (prob_state.size()>1)   // The state is stored as a vector of states/probabilities.
    {
      for (const typename LTS_TYPE::probabilistic_state_t::state_probability_pair& prob_pair : prob_state)
      {
        const state_type s = prob_pair.state();
        if (block_index_of_a_state[s] == b.key)
        {
          prob_to_block = prob_to_block + prob_pair.probability();
        }
      }
    }
    else // The state consists of a number with implicit probability 1.
    {
      const state_type s = prob_state.get();
      if (block_index_of_a_state[s] == b.key)
      {
        prob_to_block = prob_to_block + LTS_TYPE::probabilistic_state_t::probability_t::one();
      }
    }

    return prob_to_block;
  }

  typename LTS_TYPE::probabilistic_state_t calculate_new_probabilistic_state(typename LTS_TYPE::probabilistic_state_t ps)
  {
    typename LTS_TYPE::probabilistic_state_t new_prob_state;

    /* Iterate over all state probability pairs in the selected probabilistic state*/
    if (ps.size()>1)   // The state is stored as a vector of states/probabilities.
    {
      std::map<state_type, probability_fraction_type> prob_state_map;
      for (const typename LTS_TYPE::probabilistic_state_t::state_probability_pair& sp_pair : ps)
      {
        /* Check the resulting action state in the final State partition */
        state_type new_state = get_eq_class(sp_pair.state());

        if (prob_state_map.count(new_state) == 0)
        {
          /* The state is not yet in the mapping. Add the state with its probability*/
          prob_state_map[new_state] = sp_pair.probability();
        }
        else
        {
          /* The state is already in the mapping. Sum up probabilities */
          prob_state_map[new_state] = prob_state_map[new_state] + sp_pair.probability();
        }
      }
      /* Add all the state probabilities pairs in the mapping to its actual data type*/
      /* Add all the state probabilities pairs in the mapping to its actual data type*/
      typename std::map<state_type, probability_fraction_type>::const_iterator i = prob_state_map.begin();
      if (++i==prob_state_map.end()) // There is only one state with probability 1. 
      {
        new_prob_state.set(prob_state_map.begin()->first);
      }
      else
      {
        // This probabilistic state has more components. 
        for (const std::pair<const state_type, probability_fraction_type>& i: prob_state_map)
        {
          new_prob_state.add(i.first, i.second);
        }
      }

    }
    else // The state consists of a number with implicit probability 1.
    {
      /* Check the resulting action state in the final State partition */
      new_prob_state.set(get_eq_class(ps.get()));
    }

    return new_prob_state;
  }

  typename LTS_TYPE::probabilistic_state_t calculate_equivalent_probabilistic_state(step_class_type& sc)
  {
    typename LTS_TYPE::probabilistic_state_t equivalent_prob_state;

    /* Select the first probabilistic state of the step class */
    distribution_type* d = sc.distributions.front();

    typename LTS_TYPE::probabilistic_state_t old_prob_state = aut.probabilistic_state(d->key);

    equivalent_prob_state = calculate_new_probabilistic_state(old_prob_state);

    return equivalent_prob_state;
  }

  /** \brief  Two-phased partitioning algorithm described in page 204. Fig 9. Baier.
  *   \detail Refinement of state partition and step partition until no new blocks/step classes
  *           are in front of the partition lists
  */
  void refine_partition_until_it_becomes_stable()
  {
    std::list<step_class_type*> step_partition_old;
    std::list<block_type*> state_partition_old;

    // Repeat until no new blocks in front of the partition lists
    while (state_partition.front()->is_in_new_blocks == true || 
           (step_partition.size()>0 && step_partition.front()->is_in_new_step_classes == true))
    { 
      // Phase 1: Splitting of step_partition via split(M,C)
      if (state_partition.front()->is_in_new_blocks == true)
      {
        // Choose a new block in front of the state partition and change it to the back of the list
        block_type* c_block = state_partition.front();
        state_partition.pop_front();
        state_partition.push_back(c_block);
        c_block->is_in_new_blocks = false;

        // swap elements of step_partition (M) to step_partition_old (M_old)
        step_partition_old.swap(step_partition);

        // vector to add the new step classes
        static std::vector<typename std::list<step_class_type*>::iterator> pending_new_step_classes;
        pending_new_step_classes.clear();

        // iterate over all step classes
        for (typename std::list<step_class_type*>::iterator sc_iter = step_partition_old.begin(); 
          sc_iter != step_partition_old.end(); ++sc_iter)
        {
          step_class_type* sc_ptr = *sc_iter;

          // Mapping to sort the distributions based on its probability to reach a block, instead of using
          // an ordered balanced tree as suggsted in Baier
          static std::map< probability_fraction_type, std::list<distribution_type*> > distributions_ordered_by_prob;
          distributions_ordered_by_prob.clear();

          // Iterate over all distributions d of the step class and add probability to block to vector
          for (distribution_type* d : sc_ptr->distributions)
          {
            probability_fraction_type probability = probability_to_block(*d, *c_block);
            distributions_ordered_by_prob[probability].push_back(d);
          }

          // if there are multiple elements in the distributions_ordered_by_prob then we have to 
          // add them to the step partition as a new step class
          if (distributions_ordered_by_prob.size() >= 2)
          {
            step_class_type* new_step_class_ptr;

            step_classes.resize(step_classes.size() + distributions_ordered_by_prob.size() - 1);
            // iterate over mapping
            std::size_t new_class_count = 0;
            for (std::pair< probability_fraction_type, std::list<distribution_type*> > ordered_dist : distributions_ordered_by_prob)
            {
              std::list<distribution_type*>& distribution_list = ordered_dist.second;
              
              if (new_class_count == 0)
              {
                // if it is the first element of the mapping then we do not create a new step class; instead, we 
                // add its elements into the current step class that is being split.
                sc_ptr->distributions.swap(distribution_list);

                // clear the prev states of the current step class
                for (std::size_t i = 0; i < sc_ptr->prev_states.size(); i++)
                {
                  sc_ptr->prev_states[i] = false;
                }

                // recalculate prev states based on incoming transitions of each distribution
                for (distribution_type* d : sc_ptr->distributions)
                {
                  for (transition* t_ptr : d->incoming_transitions_per_label[sc_ptr->action])
                  {
                    sc_ptr->prev_states[t_ptr->from()] = true;
                  }
                }

                // add to the front of the step partition (as a new step class) if not yet added 
                if (sc_ptr->is_in_new_step_classes == false)
                {
                  // since we are iterating over this step class, we save its iterator in a pending vector,
                  // so later we can add it to the front of the step partition (outside the loop)
                  pending_new_step_classes.push_back(sc_iter);
                  sc_ptr->is_in_new_step_classes = true;
                }
              }
              else
              {
                new_step_class_ptr = &step_classes[step_classes.size() - new_class_count];

                //init new step class
                new_step_class_ptr->key = step_classes.size() - new_class_count;
                new_step_class_ptr->action = sc_ptr->action;
                new_step_class_ptr->is_in_new_step_classes = true;
                new_step_class_ptr->prev_states.resize(aut.num_states());
                new_step_class_ptr->distributions.swap(distribution_list);

                // recalculate prev states based ont incoming transitions of each distribution
                for (distribution_type* d : new_step_class_ptr->distributions)
                {
                  for (transition* t_ptr : d->incoming_transitions_per_label[new_step_class_ptr->action])
                  {
                    new_step_class_ptr->prev_states[t_ptr->from()] = true;
                  }
                }

                // add new step class to step partition and new_step_classes
                step_partition.push_front(new_step_class_ptr);
              }
              new_class_count++;
            }
          }
        }

        // Add the pending new blocks to the front of the list
        for (const typename std::list<step_class_type*>::iterator& sc_iter: pending_new_step_classes)
        {
          step_partition.splice(step_partition.begin(), step_partition_old, sc_iter);
        }

        // move remaining step classes to the end of step partition
        step_partition.splice(step_partition.end(), step_partition_old);
      }

      // Phase 2: Refinment of state_partition via Refine(X,a,M)
      if (step_partition.front()->is_in_new_step_classes == true)
      {
        // Choose some step class <a,M> in new_step_classes and remove it from new_step_classes
        step_class_type* step_class = step_partition.front();
        step_partition.pop_front();
        step_class->is_in_new_step_classes = false;
        step_partition.push_back(step_class);

        // swap elements of state_partition (X) to state_partition_old (X_old)
        state_partition_old.swap(state_partition);

        static std::vector<typename std::list<block_type*>::iterator > blocks_to_move_to_front;
        static std::vector<block_type*> new_blocks_to_move_to_front;
        blocks_to_move_to_front.clear();
        new_blocks_to_move_to_front.clear();
        // for all blocks B in X_old
        //for (block_type* b_to_split : state_partition_old)
        for (typename std::list<block_type*>::iterator block_iter = state_partition_old.begin(); 
          block_iter != state_partition_old.end(); ++block_iter)
        {
          block_type* b_to_split = *block_iter;
          block_type new_block;
          new_block.is_in_new_blocks = false;
          block_type temp_block;
          block_type* new_block_ptr;
          new_block.key = blocks.size();

          // iterate over all states in the block to split
          for (state_type s : b_to_split->states)
          {
            // if block b can reach step class <a,M> starting from state s then we move it to antoher block
            if (step_class->prev_states[s] == true)
            {
              new_block.states.push_back(s);
            }
            else
            {
              temp_block.states.push_back(s);
            }
          }

          // if both, the new block and temp block has elements, then we add a new block into the state partition
          if (new_block.states.size() > 0 && temp_block.states.size() > 0)
          {

            // First update the block_index_of_a_state by iterating over all states of the new block
            for (state_type s : new_block.states)
            {
              block_index_of_a_state[s] = new_block.key;
            }

            blocks.push_back(new_block);
            new_block_ptr = &blocks.back();

            // return states from temp_block to b_to_split
            b_to_split->states.swap(temp_block.states);
            
            // if the current block is not currently a new block then add the smaller
            // block (between new block and current block) to the list of new blocks; 
            // hence, in front of state partition
            if (b_to_split->is_in_new_blocks == false)
            {
              if (new_block_ptr->states.size() < b_to_split->states.size())
              {
                new_blocks_to_move_to_front.push_back(new_block_ptr);
                new_block_ptr->is_in_new_blocks = true;
              }
              else
              {
                b_to_split->is_in_new_blocks = true;
                blocks_to_move_to_front.push_back(block_iter);

                // add new block to the back of the state_partition list
                state_partition.push_back(new_block_ptr);
              }
            }
            else
            {
              // the current block is already in new blocks; hence, just add the new block
              // to the front of the partition as well.

              //new_blocks.push_back(new_block_ptr);
              new_block_ptr->is_in_new_blocks = true;

              // add new block to the back of the state_partition list
              new_blocks_to_move_to_front.push_back(new_block_ptr);
            }
          }
        }

        // move the remaning blocks to the end of state_partition
        state_partition.splice(state_partition.begin(), state_partition_old);

        // move the blocks that should be in new blocks to the begining of the list
        for (const typename std::list<block_type*>::iterator& block_iter: blocks_to_move_to_front)
        {
          state_partition.splice(state_partition.begin(), state_partition, block_iter);
        }

        // move the blocks that are new to the front of state_partition
        for (block_type* block_ptr : new_blocks_to_move_to_front)
        {
          state_partition.push_front(block_ptr);
        }

      }
    }
  }

  void postprocessing_stage()
  {
    //---- Post processing to keep track of the parent block of each state ----//
    //block_index_of_a_state.resize(aut.num_states());
    for (const block_type& b : blocks)
    {
      for (const state_type s : b.states)
      {
        block_index_of_a_state[s] = b.key;
      }
    }

    step_class_index_of_a_distribution.resize(aut.num_probabilistic_states());
    for (step_class_type& sc : step_classes)
    {
      sc.equivalent_step_class = sc.key;

      for (const distribution_type* d : sc.distributions)
      {
        step_class_index_of_a_distribution[d->key] = sc.key;
      }
    }

    // Merge equivalent step classes in step partition. The algorithm initially separates classes by actions,
    // but it is possible that a probabilistic sate has multiple incoming actions, hence it will be repeated
    // among multiple step classes
    typename LTS_TYPE::probabilistic_state_t new_prob_state;
    std::unordered_map<typename LTS_TYPE::probabilistic_state_t, std::vector<step_class_type*> > reduced_step_partition;

    /* Iterate over all step classes of Step_partition*/
    for (step_class_type& sc : step_classes)
    {
      if (sc.distributions.size() > 0)
      {
        typename LTS_TYPE::probabilistic_state_t new_prob_state = calculate_equivalent_probabilistic_state(sc);

        /* Add the step class index to the new probability state*/
        reduced_step_partition[new_prob_state].push_back(&sc);
      }
    }

    step_class_key_type equivalent_class_key = 0;
    for (typename std::unordered_map<typename LTS_TYPE::probabilistic_state_t,
          std::vector<step_class_type*> >::iterator i = reduced_step_partition.begin();
          i != reduced_step_partition.end(); ++i)
    {
      std::vector<step_class_type*>& sc_vector = i->second;

      equivalent_step_classes.push_back(*sc_vector[0]);
      for (step_class_type* sc :sc_vector)
      {
        sc->equivalent_step_class = equivalent_class_key;
      }
      equivalent_class_key++;
    }

  }
};


/** \brief Reduce transition system l with respect to probabilistic bisimulation.
 * \param[in/out] l The transition system that is reduced.
 * \param[in/out] timer A timer that can be used to report benchmarking results.
 */
template < class LTS_TYPE>
void probabilistic_bisimulation_reduce_bem(LTS_TYPE& l, utilities::execution_timer& timer);

/** \brief Checks whether the two initial states of two plts's are probabilistic bisimilar.
* \details This lts and the lts l2 are not usable anymore after this call.
* \param[in/out] l1 A first probabilistic transition system.
* \param[in/out] l2 A second probabilistic transition system.
* \retval True iff the initial states of the current transition system and l2 are probabilistic bisimilar */
template < class LTS_TYPE>
bool destructive_probabilistic_bisimulation_compare_bem(LTS_TYPE& l1, LTS_TYPE& l2, utilities::execution_timer& timer);


/** \brief Checks whether the two initial states of two plts's are probabilistic bisimilar.
*  \details The current transitions system and the lts l2 are first duplicated and subsequently
*           reduced modulo bisimulation. If memory space is a concern, one could consider to
*           use destructive_bisimulation_compare.
* \param[in/out] l1 A first transition system.
* \param[in/out] l2 A second transistion system.
* \retval True iff the initial states of the current transition system and l2 are probabilistic bisimilar */
template < class LTS_TYPE>
bool probabilistic_bisimulation_compare_bem(const LTS_TYPE& l1, const LTS_TYPE& l2, utilities::execution_timer& timer);


template < class LTS_TYPE>
void probabilistic_bisimulation_reduce_bem(LTS_TYPE& l, utilities::execution_timer& timer)
{

  // Apply the probabilistic bisimulation reduction algorithm.
  detail::prob_bisim_partitioner_bem<LTS_TYPE> prob_bisim_part(l,timer);

  // Clear the state labels of the LTS l
  l.clear_state_labels();

  // Assign the reduced LTS
  l.set_num_states(prob_bisim_part.num_eq_classes());
  prob_bisim_part.replace_transitions();
  prob_bisim_part.replace_probabilistic_states();
}

template < class LTS_TYPE>
bool probabilistic_bisimulation_compare_bem(
  const LTS_TYPE& l1,
  const LTS_TYPE& l2,
  utilities::execution_timer& timer)
{
  LTS_TYPE l1_copy(l1);
  LTS_TYPE l2_copy(l2);
  return destructive_probabilistic_bisimulation_compare_bem(l1_copy, l2_copy, timer);
}

template < class LTS_TYPE>
bool destructive_probabilistic_bisimulation_compare_bem(
  LTS_TYPE& l1,
  LTS_TYPE& l2,
  utilities::execution_timer& timer)
{
  std::size_t initial_probabilistic_state_key_l1;
  std::size_t initial_probabilistic_state_key_l2;

  // Merge states
  mcrl2::lts::detail::plts_merge(l1, l2);
  l2.clear(); // No use for l2 anymore.

              // The last two probabilistic states are the initial states of l2 and l1
              // in the merged plts
  initial_probabilistic_state_key_l2 = l1.num_probabilistic_states() - 1;
  initial_probabilistic_state_key_l1 = l1.num_probabilistic_states() - 2;

  detail::prob_bisim_partitioner_bem<LTS_TYPE> prob_bisim_part(l1, timer);

  return prob_bisim_part.in_same_probabilistic_class(initial_probabilistic_state_key_l2,
    initial_probabilistic_state_key_l1);
}

} // namespace mcrl2::lts::detail

#endif // MCRL2_LTS_DETAIL_LIBLTS_PBISIM_BEM_H
