// Author(s): Hector Joao Rivera Verduzco, Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file lts/detail/liblts_pbisim.h

#ifndef MCRL2_LTS_DETAIL_LIBLTS_PBISIM_GRV_H
#define MCRL2_LTS_DETAIL_LIBLTS_PBISIM_GRV_H

#include "mcrl2/lts/detail/embedded_list.h"
#include "mcrl2/lts/detail/liblts_plts_merge.h"

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/execution_timer.cxx"
#else
  import utilities;
#endif

namespace mcrl2::lts::detail
{

template < class LTS_TYPE>
class prob_bisim_partitioner_grv  // Called after Groote, Rivera Verduzco and de Vink
{
  public:
    /** \brief Creates a probabilistic bisimulation partitioner for an PLTS.
    *  \details This bisimulation partitioner applies the algorithm described in "J.F. Groote, H.J. Rivera, E.P. de Vink, 
	  *	 An O(mlogn) algorithm for probabilistic bisimulation".
    */
    prob_bisim_partitioner_grv(LTS_TYPE& l, utilities::execution_timer& timer)
      : aut(l)
    {
      mCRL2log(log::verbose) << "Probabilistic bisimulation partitioner created for " <<
                                l.num_states() << " states and " <<
                                l.num_transitions() << " transitions\n";
      timer.start("bisimulation_reduce (grv)");
      create_initial_partition();
      refine_partition_until_it_becomes_stable();
      timer.finish("bisimulation_reduce (grv)");
    }

    /** \brief Gives the number of bisimulation equivalence classes of the LTS.
    *  \return The number of bisimulation equivalence classes of the LTS.
    */
    std::size_t num_eq_classes() const
    {
      return action_constellations.size();
    }

    /** \brief Gives the bisimulation equivalence class number of a state.
    *  \param[in] A state number.
    *  \return The number of the bisimulation equivalence class to which the state belongs to. */
    std::size_t get_eq_class(const std::size_t s)
    {
      assert(s < action_states.size());
      return action_states[s].parent_block;
    }

    /** \brief Gives the bisimulation equivalence probabilistic class number of a probabilistic state.
    *  \param[in] A probabilistic state number.
    *  \return The number of the probabilistic class to which the state belongs to. */
    std::size_t get_eq_probabilistic_class(const std::size_t s)
    {
      assert(s<probabilistic_states.size());
      return probabilistic_states[s].parent_block; // The block index is the state number of the block.
    }

    /** \brief Replaces the transition relation of the current lts by the transitions
    *         of the bisimulation reduced transition system.
    * \pre The bisimulation equivalence classes have been computed. */
    void replace_transitions()
    {
      std::set<transition> resulting_transitions;

      const std::vector<transition>& trans = aut.get_transitions();
      for (const transition& t : trans)
      {
        resulting_transitions.insert(
          transition(
            get_eq_class(t.from()),
            t.label(),
            get_eq_probabilistic_class(t.to())));
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
    * \pre The bisimulation classes have been computed. */
    void replace_probabilistic_states()
    {
      std::vector<typename LTS_TYPE::probabilistic_state_t> new_probabilistic_states;

      // get the equivalent probabilistic state of each probabilistic block and add it to aut
      for (probabilistic_block_type& prob_block : probabilistic_blocks)
      {
        if (prob_block.states.size()>0)
        {
          typename LTS_TYPE::probabilistic_state_t equivalent_ps = calculate_equivalent_probabilistic_state(prob_block);
          new_probabilistic_states.push_back(equivalent_ps);
        }
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
    bool in_same_probabilistic_class_grv(const std::size_t s, const std::size_t t)
    {
      return get_eq_probabilistic_class(s) == get_eq_probabilistic_class(t);
    }


  protected:

    // --------------- BEGIN DECLARATION OF DATA TYPES ---------------------------------------------------------------

    using block_key_type = std::size_t;
    using constellation_key_type = std::size_t;
    using transition_key_type = std::size_t;
    using state_key_type = std::size_t;
    using label_type = std::size_t;
    // typedef probabilistic_arbitrary_precision_fraction probability_label_type;
    using probability_label_type = typename LTS_TYPE::probabilistic_state_t::probability_t;
    // typedef probabilistic_arbitrary_precision_fraction probability_fraction_type;
    using probability_fraction_type = typename LTS_TYPE::probabilistic_state_t::probability_t;

    struct action_transition_type : public embedded_list_node <action_transition_type>
    {
      state_key_type from;
      label_type label;
      state_key_type to;
      std::size_t* state_to_constellation_count_ptr;
    };

    struct probabilistic_transition_type : public embedded_list_node < probabilistic_transition_type >
    {
      state_key_type from = 0UL;
      probability_label_type label;
      state_key_type to = 0UL;
    };

    struct action_state_type : public embedded_list_node < action_state_type >
    {
      block_key_type parent_block = 0UL;
      std::vector<probabilistic_transition_type*> incoming_transitions;

      // Temporary
      bool mark_state = false;
      std::size_t residual_transition_cnt = 0UL;
      std::size_t* transition_count_ptr = nullptr;
    };

    struct probabilistic_state_type : public embedded_list_node < probabilistic_state_type >
    {
      block_key_type parent_block = 0UL;
      std::vector<action_transition_type*> incoming_transitions;

      // Temporary.
      bool mark_state = false;
      probability_label_type cumulative_probability;
    };

    // Prototype.
    struct action_mark_type;

    struct action_block_type : public embedded_list_node <action_block_type>
    {
      constellation_key_type parent_constellation = 0UL;
      embedded_list<action_state_type> states;
      embedded_list<probabilistic_transition_type> incoming_probabilistic_transitions;
      action_mark_type* marking;  // This value is nullptr if the block is not marked.

      action_block_type() : marking(nullptr){}
    };

    struct action_mark_type
    {
      action_block_type* action_block;
      embedded_list<action_state_type> left;
      embedded_list<action_state_type> middle;
      embedded_list<action_state_type> right;
      embedded_list<action_state_type>* large_block_ptr;

      action_mark_type(action_block_type& B) 
       : action_block(&B), 
         large_block_ptr(nullptr)
      {}
    };

    // Prototype
    struct probabilistic_mark_type;

    struct probabilistic_block_type : public embedded_list_node <probabilistic_block_type>
    {
      constellation_key_type parent_constellation = 0UL;
      embedded_list<probabilistic_state_type> states;
      probabilistic_mark_type* marking;

      // a probabilistic block has incoming action transitions ordered by label
      embedded_list<action_transition_type> incoming_action_transitions;

      probabilistic_block_type() : marking(nullptr){}
    };

    struct probabilistic_mark_type
    {
      probabilistic_block_type* probabilistic_block;
      embedded_list<probabilistic_state_type> left;
      std::vector< embedded_list<probabilistic_state_type> > middle;
      embedded_list<probabilistic_state_type> right;
      embedded_list<probabilistic_state_type>* large_block_ptr;

      probabilistic_mark_type(probabilistic_block_type& B) : probabilistic_block(&B), large_block_ptr(nullptr) {}
    };

    struct action_constellation_type 
    {
      embedded_list<action_block_type> blocks;
      std::size_t number_of_states = 0UL; // number of states in this constellation.
    }; 
    
    struct probabilistic_constellation_type 
    {
      embedded_list<probabilistic_block_type> blocks;
      std::size_t number_of_states = 0UL; // number of states in this constellation.
    };

    // --------------- END DECLARATION OF DATA TYPES ---------------------------------------------------------------

    // The class below is used to group transitions on action labels in linear space and time.
    // This makes use of the fact that action labels have a limited range, smaller than the number of transitions.
    class transitions_per_label_t
    {
      protected:
        // This is a vector that contains transitions with the same label at each entry.
        std::vector< embedded_list<action_transition_type> > m_transitions_per_label;
        // This stack indicates which positions in the vector above are occupied for efficient retrieval.
        std::stack<label_type> m_occupancy_indicator;
        
        void add_single_transition(action_transition_type& t) 
        {
          assert(t.label<m_transitions_per_label.size());
          if (m_transitions_per_label[t.label].size()==0)
          {
            m_occupancy_indicator.push(t.label);
          }
          m_transitions_per_label[t.label].push_back(t);
        }

      public:
        /* set the size of the vector m_transitions_per_label */
        void initialize(const std::size_t number_of_labels)
        { 
          m_transitions_per_label=std::vector< embedded_list<action_transition_type> >(number_of_labels);
        }

        const std::vector< embedded_list<action_transition_type> >& transitions() const
        {
          return m_transitions_per_label;
        }
        
        /* This function adds the transitions per label in order to the initial probabilistic block 
           and resets the occupancy array back to empty */
        void add_grouped_transitions_to_block(probabilistic_block_type& block)
        {
          while (!m_occupancy_indicator.empty())
          {
            const label_type action=m_occupancy_indicator.top();
            m_occupancy_indicator.pop();
            
            // The next operation resets m_transitions_per_label[action] to empty.
            assert(action<m_transitions_per_label.size());
            block.incoming_action_transitions.append(m_transitions_per_label[action]);
          } 
        }

        void move_incoming_transitions(probabilistic_state_type s, embedded_list<action_transition_type>& transition_list_with_t)
        {
          for(action_transition_type* t_ptr: s.incoming_transitions)
          {
            transition_list_with_t.erase(*t_ptr);
            add_single_transition(*t_ptr);
          }
        }

        void add_transitions(std::vector<action_transition_type>& transitions)
        {
          for(action_transition_type& t: transitions)
          {
            add_single_transition(t);
          }
        }
    };

    // The basic stores for all elementary data structures. The deque's are used intentionally
    // as there are pointers to their content. 
    std::vector<action_transition_type> action_transitions;
    std::deque<probabilistic_transition_type> probabilistic_transitions;
    std::vector<action_state_type> action_states;
    std::vector<probabilistic_state_type> probabilistic_states;
    std::deque<action_block_type> action_blocks;
    std::deque<probabilistic_block_type> probabilistic_blocks;
    std::deque<action_constellation_type> action_constellations;
    std::deque<probabilistic_constellation_type> probabilistic_constellations;

    // The storage to store the state_to_constellation counts for each transition. Transition refer with
    // a pointer to the elements in this deque. 
    std::deque<std::size_t> state_to_constellation_count;

    // The lists below contains all the non trivial constellations.
    std::stack<action_constellation_type*> non_trivial_action_constellations;
    std::stack<probabilistic_constellation_type*> non_trivial_probabilistic_constellations;
    
    // temporary data structures. declared here to prevent redeclaring these too often.
    std::deque<action_mark_type> marked_action_blocks;
    std::deque<probabilistic_mark_type> marked_probabilistic_blocks;
    std::vector< std::pair<probability_label_type, probabilistic_state_type*> > grouped_states_per_probability_in_block;
    
    // The following is a temporary data structure used to group incoming transitions on labels, which is declared
    // globally to avoid declaring it repeatedly which is time consuming.
    transitions_per_label_t transitions_per_label;


    LTS_TYPE& aut;

    bool check_data_structure()
    {
      // Check whether the constellation count in action transitions is ok. 
      for(const action_transition_type& t: action_transitions)
      {
        std::size_t constellation=probabilistic_blocks[probabilistic_states[t.to].parent_block].parent_constellation;
        std::size_t count_state_to_constellation=0;
        
        for(const action_transition_type& t1: action_transitions)
        {
          if (t.from==t1.from &&
              t.label==t1.label &&
              constellation==probabilistic_blocks[probabilistic_states[t1.to].parent_block].parent_constellation)
          {
            count_state_to_constellation++;
          }
        }
        if (count_state_to_constellation!=*t.state_to_constellation_count_ptr)
        {
          mCRL2log(log::error) << "Transition " << t.from << "--" << t.label << "->" << t.to << " has inconsistent constellation_count: " <<
                                *t.state_to_constellation_count_ptr << ". Should be " << count_state_to_constellation << ".\n";
          return false;
                                  
        }
      }

      // Check whether the number of states in an action constellation are correct.
      for(const action_constellation_type& c: action_constellations)
      {
        std::size_t counted_states=0;
        for(const action_block_type& b: c.blocks)
        {
          counted_states=counted_states+b.states.size();
        }
        assert(counted_states==c.number_of_states);  // Number of states in this action constellation does not match the number of states in its blocks.
      }

      // Check whether the number of states in a probabilistic constellation are correct.
      for(const probabilistic_constellation_type& c: probabilistic_constellations)
      {
        std::size_t counted_states=0;
        for(const probabilistic_block_type& b: c.blocks)
        {
          counted_states=counted_states+b.states.size();
        }
        assert(counted_states==c.number_of_states); // Number of states in this probabilistic constellation does not match the number of states in its blocks.
      }

      for(const action_block_type& b: action_blocks)
      {
        assert(b.states.size()>0); // Action block contains no states.

        assert(b.marking==nullptr); // Action block's marking must be a null ptr.
      }

      for(const probabilistic_block_type& b: probabilistic_blocks)
      {
        assert(b.marking==nullptr); // Probabilistic block's marking must be a null ptr.
      }

      return true;
    }

    void print_structure(const std::string& info)
    {
      std::cerr << info << " ---------------------------------------------------------------------- \n";
      std::cerr << "Number of action blocks " << action_blocks.size() << "\n";
      std::cerr << "Number of probabilistic blocks " << probabilistic_blocks.size() << "\n";
      std::cerr << "Number of action constellations " << action_constellations.size() << "\n";
      std::cerr << "Number of probabilistic constellations " << probabilistic_constellations.size() << "\n";
      for(const action_block_type b: action_blocks)
      {
        std::cerr << "ACTION BLOCK INFO ------------------------\n";
        std::cerr << "PARENT CONSTELLATION " << b.parent_constellation << "\n";
        std::cerr << "NR STATES " << b.states.size() << "\n";
        for(const probabilistic_transition_type t: b.incoming_probabilistic_transitions)
        {
          std::cerr << "INCOMING TRANS " << t.from << " --" << t.label << "-> " << t.to << "\n";
        }
      }

      for(const probabilistic_block_type b: probabilistic_blocks)
      {
        std::cerr << "probabilistic BLOCK INFO ------------------------\n";
        std::cerr << "PARENT CONSTELLATION " << b.parent_constellation << "\n";
        std::cerr << "NR STATES " << b.states.size() << "\n";
        for(const action_transition_type t: b.incoming_action_transitions)
        {
          std::cerr << "INCOMING TRANS " << t.from << " --" << t.label << "-> " << t.to << "\n";
        }
      }

    }

    /** Creates the initial partition.
    *  \details The blocks are initially partitioned based on the actions that can perform. 
    */
    void create_initial_partition()
    {
      transitions_per_label.initialize(aut.num_action_labels());
      
      // Preprocessing initialization. First we have to initialise the action/probabilistic states and
      // transitions.
      preprocessing_stage();

      // Add action transitions to its respective list grouped by label.
      transitions_per_label.add_transitions(action_transitions);

      // We start with all the action states in one block and then we refine this block
      // according to the outgoing transitions.
      action_block_type initial_action_block;

      // Link all the action states together to the initial block
      for (action_state_type& s: action_states)
      {
        s.parent_block = 0; // Initial block has number 0. 
        initial_action_block.states.push_back(s);
      }
      assert(aut.num_states()==initial_action_block.states.size());

      // Add the linked states to the list of states in the initial block
      action_blocks.push_back(initial_action_block);

      // Refine the initial action block based on the outgoing transitions.
      // This corresponds to line 5 of Algorithm 2 in [GRV].
      refine_initial_action_block(transitions_per_label.transitions());

      // Initialise the probabilistic block. Initally, there is only one block of probabilistic states.
      probabilistic_blocks.emplace_back();

      probabilistic_block_type& initial_probabilistic_block=probabilistic_blocks.back();
      initial_probabilistic_block.parent_constellation = 0;

      // Link all the probabilistic states together to the initial block
      for (probabilistic_state_type& s : probabilistic_states)
      {
        s.parent_block = 0; // The initial block has number zero.
        initial_probabilistic_block.states.push_back(s);
      }
      assert(aut.num_probabilistic_states()==initial_probabilistic_block.states.size());
      
      // Initialise the probabilistic and action constellations; they will contain
      // all the blocks.
      probabilistic_constellation_type initial_probabilistic_const;
      initial_probabilistic_const.number_of_states = aut.num_probabilistic_states();

      // Initially there is only one block in the probabilistic constellation.
      // initial_probabilistic_const.blocks.init(&probabilistic_blocks.front(), &probabilistic_blocks.front(), 1);
      initial_probabilistic_const.blocks.push_back(probabilistic_blocks.front());
      probabilistic_constellations.push_back(initial_probabilistic_const);

      // Initialise the initial action constellation.
      action_constellations.emplace_back();
      action_constellation_type& initial_action_const=action_constellations.back();

      initial_action_const.number_of_states = aut.num_states();

      // Construct the list of action blocks by linking them together.
      for(action_block_type& b : action_blocks)
      {
        b.parent_constellation = 0; // The initial constellation has number 0;
        initial_action_const.blocks.push_back(b);
      }

//-----------------------------------------------------------------------------
      // state_to_const_count_temp is used to keep track of the block to constellation count per label of each state.
      std::vector<std::size_t*> new_count_ptr(aut.num_states(),nullptr);
      
      for (const embedded_list<action_transition_type>& at_list_per_label : transitions_per_label.transitions())
      {
        // Create a new position in state_to_constellation_count if no such count exists for the parent block of the transition. 
        // Assign this position to t.state_to_constellation_count_ptr and increment its count.
        for(action_transition_type& t: at_list_per_label)
        {
          assert(t.from<new_count_ptr.size());
          if (new_count_ptr[t.from]==nullptr)
          {
            state_to_constellation_count.push_back(0);
            new_count_ptr[t.from] = &state_to_constellation_count.back();
          }
          t.state_to_constellation_count_ptr = new_count_ptr[t.from];
          (*new_count_ptr[t.from])++;
        }
    
        // Reset all the variables used to prepare to next iteration.
        for(const action_transition_type& t: at_list_per_label)
        {
          new_count_ptr[t.from] = nullptr;
        }
      } 

//-----------------------------------------------------------------------------

      // Since the transitions are already grouped by label, add them to the
      // initial probabilistic block as incoming transitions.
      transitions_per_label.add_grouped_transitions_to_block(initial_probabilistic_block);

      // Initialise the incoming probabilistic transitions for all action blocks. To that end,
      // iterate over all probabilistic transitions and add it to its respective destination block.
      for(probabilistic_transition_type& t : probabilistic_transitions)
      {
        action_state_type& s = action_states[t.to];
        action_block_type& block = action_blocks[s.parent_block];
        block.incoming_probabilistic_transitions.push_back(t);
      }

      if (initial_action_const.blocks.size()>1) 
      {
        non_trivial_action_constellations.push(&initial_action_const);
      }

      assert(check_data_structure());
    }

    /* This function performs the preprocessing stage to prepare to apply the algorithm.
       It also initialises the action and probabilistic states and transitions.  */
    void preprocessing_stage()
    {
      // Allocate space for states and transitions

      action_states.resize(aut.num_states());
      probabilistic_states.resize(aut.num_probabilistic_states());
      action_transitions.resize(aut.num_transitions());

      // Initialise the action transitions
      transition_key_type t_key = 0;
      for (const transition& t : aut.get_transitions())
      {
        action_transition_type& at = action_transitions[t_key];
        at.from = t.from();
        at.label = t.label();
        at.to = t.to();
        at.state_to_constellation_count_ptr = nullptr;

        // save incoming transition in state
        probabilistic_states[at.to].incoming_transitions.push_back(&at);

        t_key++;
      }

      // Initialise the probabilistic transitions. To this end, we have to iterate over
      // all probabilistic states.
      for (std::size_t i = 0; i < aut.num_probabilistic_states(); i++)
      {
        const typename LTS_TYPE::probabilistic_state_t& ps = aut.probabilistic_state(i);
        probabilistic_transition_type pt;
        if (ps.size()>1) // The probabilistic state is stored as a vector. 
        {
          for (const typename LTS_TYPE::probabilistic_state_t::state_probability_pair& sp_pair : ps)
          {
            pt.from = i;
            pt.label = sp_pair.probability();
            pt.to = sp_pair.state();
            probabilistic_transitions.push_back(pt);

            // save incoming transition in state
            action_states[pt.to].incoming_transitions.push_back(&probabilistic_transitions.back());
          }
        }
        else // The probabilistic state is stored as a single state number with implicit probability 1. 
        { 
          pt.from = i;
          pt.label = LTS_TYPE::probabilistic_state_t::probability_t::one();
          pt.to = ps.get();
          probabilistic_transitions.push_back(pt);
          // save incoming transition in state
          action_states[pt.to].incoming_transitions.push_back(&probabilistic_transitions.back());
        }
      }

      // End of preprocessing.
    }

    /* Refine the initial block according to its outgoing transitions.
    */
    void refine_initial_action_block(const std::vector< embedded_list<action_transition_type> >& transitions_per_label)
    {
      // Iterate over all transitions ordered by label, and refine the block.
      for (const embedded_list<action_transition_type>& t_list : transitions_per_label)
      {
        marked_action_blocks.clear();
        // The line below garbage collects marked_action_blocks to avoid covering too much memory continuously;
        static std::size_t count=0; if (count++ == 1000) { marked_action_blocks.shrink_to_fit(); count=0; }

        for(const action_transition_type& t: t_list)
        {
          action_state_type& s = action_states[t.from];
          assert(s.parent_block<action_blocks.size());
          action_block_type& parent_block = action_blocks[s.parent_block];

          // Move state s to marked states if not already added.
          if (false == s.mark_state)
          {
            // Add parent block to the list of marked blocks if not yet added
            if (parent_block.marking==nullptr) 
            {
              marked_action_blocks.emplace_back(parent_block);
              parent_block.marking=&marked_action_blocks.back();
            }

            move_list_element_back<action_state_type>(s, parent_block.states, parent_block.marking->left);
            s.mark_state = true;
          }
        }
        
        // Split the marked blocks.
        for (action_mark_type& block_marking: marked_action_blocks)
        {
          if (0 == block_marking.action_block->states.size())
          {
            // If all states in the block are marked, then return all marked states to the block.
            block_marking.action_block->states = block_marking.left;
          }
          else
          {
            // Split the block if not all states are marked.
            action_blocks.emplace_back();
            action_block_type& new_block=action_blocks.back();
            new_block.states = block_marking.left;

            // Init parent block of each state in new block.
            for(action_state_type& s: new_block.states)
            {
              s.parent_block = action_blocks.size()-1; 
            }

          }
          
          // clean mark list
          block_marking.left.clear();
          block_marking.action_block->marking=nullptr;

        }

        // Clean the marked states.
        for(const action_transition_type& t: t_list)
        {
          action_states[t.from].mark_state = false;
        }
      }

    }

    /* Move an element of a list to the back of another list.
    */
    template <typename LIST_ELEMENT>
    void move_list_element_back(LIST_ELEMENT& s, embedded_list<LIST_ELEMENT>& source_list, embedded_list<LIST_ELEMENT>& dest_list)
    {
      source_list.erase(s);
      dest_list.push_back(s);
    }

    /** \brief Move an element of a list to the front of another the list.
    *  \details
    */
    template<typename LIST_ELEMENT>
    void move_list_element_front(LIST_ELEMENT& s, embedded_list<LIST_ELEMENT>& source_list, embedded_list<LIST_ELEMENT>& dest_list)
    {
      source_list.erase(s);
      dest_list.push_front(s);
    }

    /** \brief Refine partition until it becomes stable.
    *  \details
    */
    void refine_partition_until_it_becomes_stable()
    {

      // Refine until all the constellations are trivial.
      while (!non_trivial_probabilistic_constellations.empty() || !non_trivial_action_constellations.empty())
      {
        assert(check_data_structure());

        // Refine probabilistic blocks if a non-trivial action constellation exists.
        if (!non_trivial_action_constellations.empty())
        {
          action_constellation_type* non_trivial_action_const = non_trivial_action_constellations.top();
          non_trivial_action_constellations.pop();
          assert(non_trivial_action_const->blocks.size()>=2);

// print_structure("REFINE I");
          // Choose splitter block Bc of a non-trivial constellation C, such that |Bc| <= 1/2|C|.
          // And also split constellation C into BC and C\BC in the set of constellations.
          action_block_type* Bc_ptr = choose_action_splitter(non_trivial_action_const);

          // Derive the left, right and middle sets from mark function.
          marked_probabilistic_blocks.clear();
          // The line below garbage collects marked_action_blocks to avoid covering too much memory continuously;
          static std::size_t count=0; if (count++ == 1000) { marked_probabilistic_blocks.shrink_to_fit(); count=0; }

          mark_probabilistic(*Bc_ptr, marked_probabilistic_blocks);

          // Split every marked probabilistic block based on left, middle and right.
          for (probabilistic_mark_type& B : marked_probabilistic_blocks)
          {
            // We must know whether the current constellation is already on the stack.
            bool already_on_non_trivial_constellations_stack = probabilistic_constellations[B.probabilistic_block->parent_constellation].blocks.size()>1; 

            // First return the largest of left, middle or right to the states of current processed block.
            B.probabilistic_block->states = *B.large_block_ptr;
            B.large_block_ptr->clear();

            // Split left set of the block to another block if left has states and it is not
            // the largest block.
            if (B.left.size() > 0)
            {
              split_probabilistic_block(*B.probabilistic_block, B.left);
            }

            // Split right set of the block to another block if right has states and it is not
            // the largest block.
            if (B.right.size() > 0)
            {
              split_probabilistic_block(*B.probabilistic_block, B.right);
            }

            // Iterate over all middle sets. Split middle sets of the current block to another block if 
            //  the middle set has states and it is not the largest block.
            for (embedded_list<probabilistic_state_type>& middle : B.middle)
            {
              if (middle.size() > 0)
              {
                split_probabilistic_block(*B.probabilistic_block, middle);
              }
            }

            // Move the parent constellation of the current block to the front of the
            // constellation list if it became unstable.
            if (!already_on_non_trivial_constellations_stack &&  probabilistic_constellations[B.probabilistic_block->parent_constellation].blocks.size()>1)
            {
              probabilistic_constellation_type& parent_const = probabilistic_constellations[B.probabilistic_block->parent_constellation];
              assert(parent_const.blocks.size()>1);
              non_trivial_probabilistic_constellations.push(&parent_const);
            }

            // Reset middle vector to prepare for the next mark process
            B.middle.clear();
          }
        }

        // Refine action blocks if a non-trivial probabilistic constellation exists.
        if (!non_trivial_probabilistic_constellations.empty())
        {
// print_structure("REFINE II");

          probabilistic_constellation_type* non_trivial_probabilistic_const = non_trivial_probabilistic_constellations.top();
          non_trivial_probabilistic_constellations.pop();
          assert(non_trivial_probabilistic_const->blocks.size()>=2);
          // Choose splitter block Bc of a non-trivial constellation C, such that |Bc| <= 1/2|C|.
          // And also split constellation C into BC and C\BC in the set of constellations.
          probabilistic_block_type* Bc_ptr = choose_probabilistic_splitter(non_trivial_probabilistic_const);
          
          // For all incoming labeled "a" transitions of each state in BC call the mark function and split the blocks.
          for (typename embedded_list<action_transition_type>::iterator i=Bc_ptr->incoming_action_transitions.begin(); 
                        i!=Bc_ptr->incoming_action_transitions.end() ;  )
          {
            // Derive the left, right and middle sets from mark function based on the incoming labeled "a" transitions.
            const label_type a = i->label;
            marked_action_blocks.clear();
            mark_action(marked_action_blocks, a, i, Bc_ptr->incoming_action_transitions.end());  // The iterator i is implicitly increased 
                                                                                          // to the position in the list with the next action.

            // Split every marked probabilistic block based on left, middle and right.
            for (action_mark_type& B : marked_action_blocks)
            {
              // We must know whether the current constellation is already on the stack.
              bool already_on_non_trivial_constellations_stack = action_constellations[B.action_block->parent_constellation].blocks.size()>1; 

              // First return the largest of left, middle or right to the states of current processed block.
              B.action_block->states = *B.large_block_ptr;
              B.large_block_ptr->clear();

              // Split left set of the block to another block if left has states and it is not the largest block.
              if (B.left.size() > 0)
              {
                split_action_block(*B.action_block, B.left);
              }

              // Split right set of the block to another block if right has states and it is not the largest block.
              if (B.right.size() > 0)
              {
                split_action_block(*B.action_block, B.right);
              }

              // Split middle set of the block to another block if middle has states and it is not
              // the largest block.
              if (B.middle.size() > 0)
              {
                split_action_block(*B.action_block, B.middle);
              }

              // Move the parent constellation of the current block to the front of the
              // constellation list if it became unstable.
              if (!already_on_non_trivial_constellations_stack && action_constellations[B.action_block->parent_constellation].blocks.size()>1)
              {
                action_constellation_type& parent_const = action_constellations[B.action_block->parent_constellation];
                assert(parent_const.blocks.size()>1);
                non_trivial_action_constellations.push(&parent_const);
              }

            }
          }
        }
      }
      // print_structure("END REFINE");
      assert(check_data_structure());
    }

    /** \brief Split a probabilistic block.
    *  \details Creates another block containing the states specified in states_of_new_block. It adds the new block
    *           to the same constellation as the current block to split.
    */
    void split_probabilistic_block(probabilistic_block_type& block_to_split, embedded_list<probabilistic_state_type>& states_of_new_block)
    {
      // First create the new block to be allocated, and initialise its parameters
      probabilistic_blocks.emplace_back();
      probabilistic_block_type& new_block=probabilistic_blocks.back();
      
      new_block.parent_constellation = block_to_split.parent_constellation; // The new block is in the same constellation as B
      new_block.states = states_of_new_block;
      states_of_new_block.clear();

      // Add the incoming action transition of the new block. To this end, iterate over all
      // states in the new block and add the incoming transitions of each state to the
      // incoming transitions of the new block. Also keep track of the labels of the incoming
      // transitions.
      for(probabilistic_state_type& s: new_block.states)
      {
        // Update the parent block of the state
        s.parent_block = probabilistic_blocks.size()-1; 
        transitions_per_label.move_incoming_transitions(s,block_to_split.incoming_action_transitions);  
      }

      // Since the transitions are already grouped by label, add them to the
      // initial probabilistic block as incoming transitions.
      transitions_per_label.add_grouped_transitions_to_block(new_block);

      // Add the new block to the back of the list of blocks in the parent constellation.
      probabilistic_constellation_type& parent_const = probabilistic_constellations[new_block.parent_constellation];

      parent_const.blocks.push_back(probabilistic_blocks.back());
    }

    /** \brief Split an action block.
    *  \details Creates another block containing the states specified in states_of_new_block. It adds the new block
    *           to the same constellation as the current block to split.
    */
    void split_action_block(action_block_type& block_to_split, embedded_list<action_state_type>& states_of_new_block)
    {
      // First create the new block to be allocated, and initialise its parameters
      action_blocks.emplace_back();
      action_block_type& new_block=action_blocks.back();

      new_block.parent_constellation = block_to_split.parent_constellation; // The new block is in the same constellation as block to split
      new_block.states = states_of_new_block;
      states_of_new_block.clear();

      // Add the incoming action transition of the new block. To this end, iterate over all
      // states in the new block and add the incoming transitions of each state to the
      // incoming transitions of the new block.
      for(action_state_type& s: new_block.states)
      {
        // Update the parent block of the state
        s.parent_block = action_blocks.size()-1;   

        // Iterate over all incoming transitions of the state, to add them to the new block
        for (probabilistic_transition_type* t : s.incoming_transitions)
        {
          // Move transition from list of transitions of previous block to new block
          move_list_element_back((*t), block_to_split.incoming_probabilistic_transitions,
            new_block.incoming_probabilistic_transitions);
        }
      }

      // Add the new block to the back of the list of blocks in the parent constellation.
      action_constellation_type& parent_const = action_constellations[new_block.parent_constellation];

      parent_const.blocks.push_back(action_blocks.back());
    }

    /** \brief Gives the probabilistic blocks that are marked by block Bc.
    *  \details Derives the left, middle and rigth sets of the marked probabilistic blocks, based on the
    *           incoming probabilistic transitions in block Bc.
    */
    void mark_probabilistic(const action_block_type& Bc, std::deque<probabilistic_mark_type>& marked_probabilistic_blocks)
    {
      // First, iterate over all incoming transitions of block Bc. Mark the blocks that are reached
      // and calculate the cumulative probability of the reached state. This is the probability of
      // a state to reach block Bc.
      for( probabilistic_transition_type& pt: Bc.incoming_probabilistic_transitions)
      {
        probabilistic_state_type& s = probabilistic_states[pt.from];
        const probability_label_type& p = pt.label;
        probabilistic_block_type& B = probabilistic_blocks[s.parent_block];

        // If the block was not previously marked, then mark the block and move all states to right
        if (nullptr == B.marking)
        {
          marked_probabilistic_blocks.emplace_back(B);
          B.marking = &marked_probabilistic_blocks.back();
          B.marking->right = B.states;
          B.states.clear(); 

          // Also initialise the larger block pointer to the right set and the maximum cumulative
          // probability of the block to p.
          B.marking->large_block_ptr = &B.marking->right;
        }

        // Since state s can reach block Bc, move state s to the left set if not yet added, and mark the state
        if (false == s.mark_state)
        {
          // State s is not yet marked. Mark the state and move it to left set. In addition, initialise
          // its cumulative probability.
          s.mark_state = true;
          s.cumulative_probability = p;
          move_list_element_back<probabilistic_state_type>(s, B.marking->right, B.marking->left);
        }
        else 
        {
          // State s was already added to left. Just update its cumulative probability.
          s.cumulative_probability = s.cumulative_probability + p;
        }
      }

      // Group all states with the same cumulative probability to construct the middle sets.
      // To this end, iterate over all marked blocks. For each block, first add all the states
      // with probability lower than the max_cumulative_probability of the block to the middle set.
      for (probabilistic_mark_type& B : marked_probabilistic_blocks)
      {
        // Clear this locally used data structure and reset it so now and then to avoid that it requires
        // too much data. 
        grouped_states_per_probability_in_block.clear();
        static std::size_t count=0; if (count++ == 1000) { marked_action_blocks.shrink_to_fit(); count=0; }

        embedded_list<probabilistic_state_type> middle_temp=B.left;
        B.left.clear();

        // First, add all states lower than max_cumulative_probability to the middle set.
        for(typename embedded_list<probabilistic_state_type>::iterator i=middle_temp.begin(); i!=middle_temp.end(); )
        {
          probabilistic_state_type& s= *i;
          i++; // Increment the iterator here, such that we can change the list. 

          if (s.cumulative_probability == probability_label_type().one()) 
          {
            // State s has probability lower than max_cumulative_probability. Add to the middle set.
            move_list_element_back<probabilistic_state_type>((s), middle_temp, B.left);
          }

          // Also reset the marked_state variable in the state here, taiking advantage that we
          // are iterating over all marked states
          s.mark_state = false;
        }

        // Add all the states corresponding to the bigger and smaller probability to middle.
        // Save the remaining states in a vector to sort them later.
        for(probabilistic_state_type& s: middle_temp)
        {
          grouped_states_per_probability_in_block.emplace_back(s.cumulative_probability, &s);
        }

        // Sort the probabilities of middle, not including the biggest and smallest probability
        std::sort(grouped_states_per_probability_in_block.begin(), grouped_states_per_probability_in_block.end());
        
        // Construct the rest of the middle set based on the grouped probabilities. To this end, 
        // traverse all the vector with the grouped states by probability. Store the states with
        // the same probability to a new sub-set in middle set. current_probability is used to
        // keep track of the probability that is currently being processed.
        probability_label_type current_probability = probability_label_type().zero();

        if (grouped_states_per_probability_in_block.size()>0)
        { 
          B.middle.emplace_back();
          for (const std::pair<probability_label_type, probabilistic_state_type*>& cumulative_prob_state_pair : grouped_states_per_probability_in_block)
          {
            probabilistic_state_type* s = cumulative_prob_state_pair.second;
            if (current_probability != cumulative_prob_state_pair.first)
            {
              // The current state has a different probability as the current probability. Allocate
              // another space in the middle to store this state and change the current probability.
              current_probability = cumulative_prob_state_pair.first;
              B.middle.emplace_back(); //put empty list at end of B.middle.
            }

            move_list_element_back<probabilistic_state_type>((*s), middle_temp, B.middle.back());
          }
        }

        // Now that we have all the states in left, middle and right; we have to find the largest
        // set. The large block is initialised to be the right set; hence, we only compare if
        // left and middle are larger.
        if (B.left.size() > B.large_block_ptr->size())
        {
          B.large_block_ptr = &B.left;
        }

        // Iterate over all subsets of middle set to see if there is a one that is the largest.
        for (embedded_list<probabilistic_state_type>& middle_set : B.middle)
        {
          if (middle_set.size() > B.large_block_ptr->size())
          {
            B.large_block_ptr = &middle_set;
          }
        }

        // Finally, reset the block_is_marked variable of the current block.
        B.probabilistic_block->marking = nullptr; 
      }
    }

    /** \brief Gives the action blocks that are marked by probabilistic block Bc.
    *  \details Derives the left, middle and rigth sets of the marked action blocks, based on the
    *           incoming action transitions labeled with "a" in block Bc.
    */
    void mark_action(std::deque<action_mark_type>& marked_action_blocks, 
                     const label_type& a, 
                     typename embedded_list<action_transition_type>::iterator& action_walker_begin,
                     const typename embedded_list<action_transition_type>::iterator action_walker_end)
    {
      assert(action_walker_begin!=action_walker_end && action_walker_begin->label==a);

      // For all incoming transitions with label "a" of block Bc calculate left, middle and right.
      // To this end, first move all the states of the block that was reached by traversing the
      // transition backwards to its right set, then move all the states that can reach block Bc with 
      // an "a" action to left and decrement the residual transition count of the state.
      for(typename embedded_list<action_transition_type>::iterator action_walker=action_walker_begin; 
          action_walker!=action_walker_end && action_walker->label==a;  
          action_walker++)
      {
        action_transition_type& t= *action_walker;
        action_state_type& s = action_states[t.from];
        action_block_type& B = action_blocks[s.parent_block];  
        
        // If the block was not previously marked, then mark the block and add all states to right.
        if (nullptr == B.marking)
        {
          marked_action_blocks.emplace_back(B);
          B.marking = &marked_action_blocks.back();
          B.marking->right=B.states;
          B.states.clear();
          // Also initialise the larger block pointer to the right set.
          B.marking->large_block_ptr = &B.marking->right;
        }

        // Since state s can reach block Bc, move state s to the left set if not yet added, and mark the state
        if (false == s.mark_state)
        {
          // State s is not yet marked. Mark the state and move it to left set. In addition, initialise
          // its residual transition count.
          s.mark_state = true;
          s.residual_transition_cnt = *t.state_to_constellation_count_ptr;
          move_list_element_back<action_state_type>(s, B.marking->right, B.marking->left);
        }

        s.residual_transition_cnt--;
      }

      // Now, for all marked blocks, check the residual transition count of all the states in left. If the transition 
      // count is zero, it means that the state only can reach block BC. If the transition count is greater than 
      // zero, the state has transitions to the other part of the constellation; hence, those states have to be
      // moved to middle.
      for (action_mark_type& B : marked_action_blocks)
      {
        // Iterate over all left states in B and check whether the state has to be moved to middle.
        for(typename embedded_list<action_state_type>::iterator i=B.left.begin(); i!=B.left.end(); )
        {
          action_state_type& s= *i; 
          i++;  // This iterator is incremented here as s will be removed from the list over which iteration takes place.
          if (s.residual_transition_cnt > 0)
          {
            // The transition count is greater than zero. Move state to middle set.
            move_list_element_back<action_state_type>(s, B.left, B.middle);
          }

          // Also reset the marked_state variable in the state here, taking advantage that we
          // are iterating over all marked states
          s.mark_state = false;
        }

        // Find the largest set.
        if (B.left.size() > B.large_block_ptr->size())
        {
          B.large_block_ptr = &B.left;
        }
        if (B.middle.size() > B.large_block_ptr->size())
        {
          B.large_block_ptr = &B.middle;
        }

        // Finally, reset the block_is_marked variable of the current block.
        B.action_block->marking = nullptr;
      }

      // Update the state_to_constellation_count of each transition. Increment action_walker_begin
      // such that it points to the next action after this loop. 
      for( ;
          action_walker_begin!=action_walker_end && action_walker_begin->label==a;  
          action_walker_begin++)
      {
        action_transition_type& t= *action_walker_begin;
        action_state_type& s = action_states[t.from];

        // If the residual_transition_cnt is greater than zero, it means that the state
        // is in the middle set; hence, the state_to_constellation_count has to be updated.
        if (s.residual_transition_cnt > 0)
        {
          std::size_t state_to_constellation_count_old = *t.state_to_constellation_count_ptr;

          if (state_to_constellation_count_old != s.residual_transition_cnt)
          {
            // This is the first transition from this state to a new block. 
            // First update the state_to_constellation_count with the residual_transition_cnt 
            // which is used by the transitions that we do not visit. Also the not yet 
            // visited transitions are set to this value.
            *t.state_to_constellation_count_ptr = s.residual_transition_cnt;

            // Now allocate another state_to_constellation_count for the Bc block
            state_to_constellation_count.emplace_back(state_to_constellation_count_old - s.residual_transition_cnt);
            s.transition_count_ptr = &state_to_constellation_count.back();
          }
          t.state_to_constellation_count_ptr = s.transition_count_ptr; 
        }
      }
    }

    /** \brief Choose an splitter block from a non trivial constellation.
    *  \details The number of states in the chosen splitter is always less than the half of the non 
    *           trivial constellation. Furtheremore, the selected splitter is moved to a new constellation.
    */
    action_block_type* choose_action_splitter(action_constellation_type* non_trivial_action_const)
    {
      assert(non_trivial_action_const->blocks.size()>=2);

      // First, determine the block to split from constellation. It is the block |Bc| < 1/2|C|
      // First try with the first block in the list.
      action_block_type* Bc = &non_trivial_action_const->blocks.front();

      if (Bc->states.size() > (non_trivial_action_const->number_of_states / 2))
      {
        // The block is bigger than 1/2|C|. Choose another one.
        Bc = &non_trivial_action_const->blocks.back();
      }

      // Now split the block of the constellation.
      // First unlink Bc from the list of blocks of the non trivial constellation.
      non_trivial_action_const->blocks.erase(*Bc);

      // Update the number of states and blocks of the non trivial block
      non_trivial_action_const->number_of_states -= Bc->states.size();

      // Check if the constellation is still non-trivial; if not, move it to the non trivial constellation stack.
      if (non_trivial_action_const->blocks.size() > 1)  
      {
        //The constellation is non trivial, put it in the stack of non_trivial_constellations.
        non_trivial_action_constellations.push(non_trivial_action_const);
      } 

      // Add Bc to a new constellation
      action_constellations.emplace_back();
      action_constellation_type& new_action_const=action_constellations.back();

      Bc->parent_constellation = action_constellations.size()-1;
      new_action_const.blocks.push_back(*Bc);
      new_action_const.number_of_states = Bc->states.size();

      return Bc;
    }

    /** \brief Choose an splitter block from a non trivial constellation.
    *  \details The number of states in the chosen splitter is always less than the half of the non
    *           trivial constellation. Furtheremore, the selected splitter is moved to a new constellation.
    */
    probabilistic_block_type* choose_probabilistic_splitter(probabilistic_constellation_type* non_trivial_probabilistic_const)
    {
      // First, determine the block to split from constellation. It is the block |Bc| < 1/2|C|
      // First try with the first block in the list.
      probabilistic_block_type* Bc = &non_trivial_probabilistic_const->blocks.front();
      if (Bc->states.size() > (non_trivial_probabilistic_const->number_of_states / 2))
      {
        // The block is bigger than 1/2|C|. Choose another one.
        Bc = &non_trivial_probabilistic_const->blocks.back();
      }

      // Now split the block of the constellation.
      // First unlink Bc from the list of blocks of the non trivial constellation.
      non_trivial_probabilistic_const->blocks.erase(*Bc);

      // Update the number of states and blocks of the non trivial block
      non_trivial_probabilistic_const->number_of_states -= Bc->states.size();

      // Check if the constellation is still non-trivial; if not, move it to the trivial constellation stack.
      if (non_trivial_probabilistic_const->blocks.size() > 1)
      {
        //The constellation is non trivial, put it in the stack of non_trivial_constellations.
        non_trivial_probabilistic_constellations.push(non_trivial_probabilistic_const);
      }

      // Add Bc to a new constellation
      probabilistic_constellations.emplace_back();
      probabilistic_constellation_type& new_probabilistic_const=probabilistic_constellations.back();
      Bc->parent_constellation = probabilistic_constellations.size()-1;
      new_probabilistic_const.blocks.push_back(*Bc); 
      new_probabilistic_const.number_of_states = Bc->states.size();

      return Bc;
    }

    typename LTS_TYPE::probabilistic_state_t calculate_new_probabilistic_state(typename LTS_TYPE::probabilistic_state_t ps)
    {
      typename LTS_TYPE::probabilistic_state_t new_prob_state;

      /* Iterate over all state probability pairs in the selected probabilistic state*/
      if (ps.size()>1)  // The state is stored as a vector of states and probabilities.
      {
        std::map <state_key_type, probability_fraction_type> prob_state_map;
        for (const typename LTS_TYPE::probabilistic_state_t::state_probability_pair& sp_pair : ps)
        {
          /* Check the resulting action state in the final State partition */
          state_key_type new_state = get_eq_class(sp_pair.state());

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
        typename std::map<state_key_type, probability_fraction_type>::iterator i = prob_state_map.begin();
        if (++i==prob_state_map.end()) // There is only one state with probability one. 
        {
          new_prob_state.set(prob_state_map.begin()->first);
        }
        else
        {
          // This probabilistic state has more components. 
          for (const auto& i: prob_state_map)
          {
            new_prob_state.add(i.first, i.second);
          }
        }
      }
      else // The state is stored as a number with implicit probability 1. 
      {
        /* Check the resulting action state in the final State partition */
        new_prob_state.set(get_eq_class(ps.get()));
      }
      return new_prob_state;
    }

    typename LTS_TYPE::probabilistic_state_t calculate_equivalent_probabilistic_state(probabilistic_block_type& pb)
    {
      // Select the first probabilistic state of the probabilistic block.
      assert(pb.states.size()>0);
      const probabilistic_state_type& s = pb.states.front();

      if (s.incoming_transitions.size()>0)
      {
        // Take an incoming transition to know the key of the state
        state_key_type s_key = s.incoming_transitions.back()->to;

        const typename LTS_TYPE::probabilistic_state_t& old_prob_state = aut.probabilistic_state(s_key);

        return calculate_new_probabilistic_state(old_prob_state);
      }
      // else it is the initial proabilistic state.
      return calculate_new_probabilistic_state(aut.initial_probabilistic_state());

    }
};


/** \brief Reduce transition system l with respect to probabilistic bisimulation.
* \param[in/out] l The transition system that is reduced.
*/
template < class LTS_TYPE>
void probabilistic_bisimulation_reduce_grv(LTS_TYPE& l, utilities::execution_timer& timer);

/** \brief Checks whether the two initial states of two plts's are probabilistic bisimilar.
* \details This lts and the lts l2 are not usable anymore after this call.
* \param[in/out] l1 A first probabilistic transition system.
* \param[in/out] l2 A second probabilistic transition system.
* \retval True iff the initial states of the current transition system and l2 are probabilistic bisimilar */
template < class LTS_TYPE>
bool destructive_probabilistic_bisimulation_compare_grv(LTS_TYPE& l1, LTS_TYPE& l2, utilities::execution_timer& timer);

/** \brief Checks whether the two initial states of two plts's are probabilistic bisimilar.
*  \details The current transitions system and the lts l2 are first duplicated and subsequently
*           reduced modulo bisimulation. If memory space is a concern, one could consider to
*           use destructive_bisimulation_compare.
* \param[in/out] l1 A first transition system.
* \param[in/out] l2 A second transistion system.
* \retval True iff the initial states of the current transition system and l2 are probabilistic bisimilar */
template < class LTS_TYPE>
bool probabilistic_bisimulation_compare_grv(const LTS_TYPE& l1, const LTS_TYPE& l2, utilities::execution_timer& timer);

template < class LTS_TYPE>
void probabilistic_bisimulation_reduce_grv(LTS_TYPE& l, utilities::execution_timer& timer)
{
  // Apply the probabilistic bisimulation reduction algorithm.
  detail::prob_bisim_partitioner_grv<LTS_TYPE> prob_bisim_part(l, timer);

  // Clear the state labels of the LTS l
  l.clear_state_labels();

  // Assign the reduced LTS
  l.set_num_states(prob_bisim_part.num_eq_classes());
  prob_bisim_part.replace_transitions();
  prob_bisim_part.replace_probabilistic_states();
}

template < class LTS_TYPE>
bool probabilistic_bisimulation_compare_grv(
  const LTS_TYPE& l1,
  const LTS_TYPE& l2,
  utilities::execution_timer& timer)
{
  LTS_TYPE l1_copy(l1);
  LTS_TYPE l2_copy(l2);
  return destructive_probabilistic_bisimulation_compare_grv(l1_copy, l2_copy, timer);
}

template < class LTS_TYPE>
bool destructive_probabilistic_bisimulation_compare_grv(
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
  // in the merged plts.
  initial_probabilistic_state_key_l2 = l1.num_probabilistic_states() - 1;
  initial_probabilistic_state_key_l1 = l1.num_probabilistic_states() - 2;

  detail::prob_bisim_partitioner_grv<LTS_TYPE> prob_bisim_part(l1,timer);

  return prob_bisim_part.in_same_probabilistic_class_grv(initial_probabilistic_state_key_l2, 
                                                         initial_probabilistic_state_key_l1);
}

} // end namespace detail
// end namespace lts
// end namespace mcrl2
#endif //MCRL2_LTS_DETAIL_LIBLTS_PBISIM_GRV_H
