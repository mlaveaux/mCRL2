// Author(s): Muck van Weerdenburg
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file sim.cpp

#include "mcrl2/utilities/input_tool.h"
#include "mcrl2/data/rewriter_tool.h"
#include "mcrl2/lps/io.h"
#include "mcrl2/lps/detail/instantiate_global_variables.h"
#include "mcrl2/lts/simulation.h"

using namespace mcrl2;
using namespace mcrl2::utilities;
using namespace mcrl2::utilities::tools;
using namespace mcrl2::core;
using namespace mcrl2::core::detail;
using namespace mcrl2::data;
using namespace mcrl2::lps;

using mcrl2::data::tools::rewriter_tool;

class sim_tool : public rewriter_tool< input_tool >
{
    typedef rewriter_tool<input_tool> super;

  private:
    static std::string print_state(const state& s)
    {
      std::string output;
      bool first = true;
      for (const data_expression& d: s)
      {
        if (!first)
        {
          output += ", ";
        }
        first = false;

        if (mcrl2::data::is_variable(d))
        {
          output += "_";
        }
        else
        {
          output += data::pp(d);
        }   
      }       
      return "[" + output + "]";
    }     

    static std::string print_state(const stochastic_state& s)
    {
      if (s.size()>1)
      {
        return "[ probabilistic state of size " + std::to_string(s.size()) + " ]";
      }

      return print_state(s.states[0]);
    }

  protected:

    bool m_do_not_use_dummies;

    void add_options(interface_description& desc)
    {
      super::add_options(desc);
      desc.add_option("nodummy", "do not replace global variables in the LPS with dummy values", 'y');
    }

    void parse_options(const command_line_parser& parser)
    {
      super::parse_options(parser);
      m_do_not_use_dummies = 0 < parser.options.count("nodummy");
    }

    void show_state_or_transition_menu(
                       const lps::simulation& simulation,
                       const std::size_t state_index)
    {
      assert(state_index < simulation.trace().size());
      const simulation::simulator_state_t& current_state = simulation.trace()[state_index];

      if (simulation.probabilistic_state_must_be_selected(state_index))
      {
        std::cout << std::endl << (state_index == 0 ? "initial" : "current") 
                  << " state: " << std::endl;

        assert(current_state.source_state.size()>0);
        for (std::size_t i = 0; i < current_state.source_state.size(); i++)
        {
          std::cout << i << ": " 
                    << ((current_state.source_state.probabilities[i]!=data::sort_real::real_one())?
                               lps::pp(current_state.source_state.probabilities[i]) + ":  ":"")
                    << print_state(current_state.source_state.states[i]) << std::endl;
        }
      }
      else // a normal state must be selected. 
      {
        std::cout << std::endl << (state_index == 0 ? "initial" : "current") 
                  << " state: " << print_state(current_state.source_state) << std::endl;

        for (std::size_t i = 0; i < current_state.transitions.size(); i++)
        {
          std::cout << i << ": " << lps::pp(current_state.transitions[i].action) << "  ->  " 
                    << print_state(current_state.transitions[i].state) << std::endl;
        }

        if (current_state.transitions.size() == 0)
        {
          std::cout << "deadlock" << std::endl;
        }
      }
    }

    void handle_help() const
    {
      std::cout <<
                "During the simulation the following commands are accepted (short/long):\n"
                "   N                take transition N to the corresponding state\n"
                "                    (where N is a number)\n"
                "   u/undo           go to previous state in trace\n"
                "   r/redo           go to next state in trace\n"
                "   i/initial        go to initial state (preserving trace)\n"
                "   g/goto N         go to position N in trace\n"
                "   t/trace          print trace (current state is indicated with '>')\n"
                "   l/load FILENAME  load trace from file FILENAME\n"
                "   s/save FILENAME  save trace to file FILENAME\n"
                "   h/help           print this help message\n"
                "   q/quit           quit\n";
    }

    void handle_number(const std::string& s,
                       lps::simulation& simulation,
                       std::size_t& state_index)
    {
      std::size_t index;
      std::stringstream ss(s);
      ss >> index;
      const simulation::simulator_state_t& current_state = simulation.trace()[state_index];

      if (simulation.probabilistic_state_must_be_selected(state_index)) // Select a probabilistic state
      {
        if (index < current_state.source_state.size())
        {
          simulation.truncate(state_index, false);
          simulation.select_state(index);
          std::cout << std::endl << "selected state: [" << lps::pp(current_state.source_state.states[index]) << "]" << std::endl;
          show_state_or_transition_menu(simulation, state_index);
        }
        else
        {
          std::cout << "invalid state index";
          if (current_state.source_state.size() > 0)
          {
            std::cout << " " << index << " (maximum is " << current_state.source_state.size() - 1 << ")";
          }
          std::cout << std::endl;
        }
      }
      else // Select a transition.
      {
        if (index < current_state.transitions.size())
        {
          simulation.truncate(state_index, false); 
          simulation.select_transition(index);
          std::cout << std::endl << "selected transition: " << lps::pp(current_state.transitions[index].action) << std::endl;
          state_index++;
          show_state_or_transition_menu(simulation, state_index);
        }
        else
        {
          std::cout << "invalid transition index";
          if (current_state.transitions.size() > 0)
          {
            std::cout << " " << index << " (maximum is " << current_state.transitions.size() - 1 << ")";
          }
          std::cout << std::endl;
        }
      }
    }

    void handle_undo(lps::simulation& simulation,
                     std::size_t& state_index)
    {
      if (state_index > 0)
      {
        state_index--;
        show_state_or_transition_menu(simulation, state_index);
      }
      else
      {
        std::cout << "already at start of trace" << std::endl;
      }
    }

    void handle_redo(lps::simulation& simulation,
                     std::size_t& state_index)
    {    
      if (state_index + 1 < simulation.trace().size())
      {
        const simulation::simulator_state_t& current_state = simulation.trace()[state_index];
        std::cout << std::endl << "transition: " << lps::pp(current_state.transitions[current_state.transition_number].action) << std::endl;
        state_index++;
        show_state_or_transition_menu(simulation, state_index);
      }
      else
      {
        std::cout << "already at end of trace" << std::endl;
      }
    }

    void handle_goto(const std::string& s,
                     lps::simulation& simulation,
                     std::size_t& state_index)
    {    
      std::istringstream sin(((s[1] == ' ') ? s.substr(2) : s.substr(5)));
      std::size_t index;
      sin >> index;
      if (index < simulation.trace().size())
      {
        state_index = index;
        show_state_or_transition_menu(simulation, state_index);
      }
      else
      {
        std::cout << "invalid trace position " << index << " (maximum is " 
                  << (simulation.trace()[state_index].transitions.size() - 1) << ")" << std::endl;
      }
    }

    void handle_trace(const lps::simulation& simulation,
                      const std::size_t& state_index)
    {    
      std::cout << std::endl << "current trace:" << std::endl;
      for (std::size_t i = 0; i < simulation.trace().size(); i++)
      {
        if (i == state_index)
        {
          std::cout << ">";
        }
        else
        {
          std::cout << " ";
        }

        if (i == 0)
        {
          std::cout << "initial state: ";
        }
        else
        {
          std::cout << lps::pp(simulation.trace()[i - 1].transitions[simulation.trace()[i - 1].transition_number].action) << "  ->  ";
        }

        std::cout << print_state(simulation.trace()[i].source_state) << std::endl;
      }
      show_state_or_transition_menu(simulation, state_index);
    }

    void handle_save(const std::string& s,
                     const lps::simulation& simulation)
    {
      std::string filename((s[1] == ' ') ? s.substr(2) : s.substr(5));
      try
      {
        simulation.save(filename);
        std::cout << "trace saved" << std::endl;
      }
      catch (mcrl2::runtime_error& err)
      {
        std::cout << "error saving trace: " << err.what() << std::endl;
      }
    }

    void handle_load(const std::string& s,
                     lps::simulation& simulation,
                     std::size_t& state_index)
    {
      std::string filename((s[1] == ' ') ? s.substr(2) : s.substr(5));
      try
      {
        simulation.load(filename);
        state_index = 0;
        std::cout << "trace loaded" << std::endl;
        show_state_or_transition_menu(simulation, state_index);
      }
      catch (mcrl2::runtime_error& err)
      {
        std::cout << "error loading trace: " << err.what() << std::endl;
      }
    }

  public:

    sim_tool()
      : super(
        "lpssim",
        "Muck van Weerdenburg",
        "Command-line simulation of an LPS",
        "Simulate the LPS in INFILE via a text-based interface."
      ),
      m_do_not_use_dummies(false)
    {}

    bool run()
    {
      mcrl2::lps::stochastic_specification spec;
      load_lps(spec, input_filename());

      if (!m_do_not_use_dummies)
      {
        lps::detail::instantiate_global_variables(spec);
      }


      // Print the parameters of the linear process.
      std::cout << "process parameters: ";
      bool first = true;
      for (const data::data_expression& d: spec.process().process_parameters())
      {
        if (!first)
        {
          std::cout << ", ";
        }

        std::cout << d << ":" << d.sort();
        first = false;
      }

      lps::simulation simulation(spec, rewrite_strategy());
      std::size_t state_index = 0;


      bool notdone = true;
      show_state_or_transition_menu(simulation, state_index);
      while (notdone)
      {
        std::string s;

        std::cout << "? ";
        std::cout.flush();
        getline(std::cin, s);
        if ((s.length() > 0) && (s[s.length()-1] == '\r'))
        {
          s.resize(s.length()-1);
        }

        if (std::cin.eof() || (s == "q") || (s == "quit"))
        {
          if (std::cin.eof())
          {
            std::cout << std::endl;
          }
          notdone = false;
        }
        else if ((s == "h") || (s == "help"))
        {
          handle_help();
        }
        else if (isdigit(s[0]))
        {
          handle_number(s, simulation, state_index);
        }
        else if ((s == "i") || (s == "initial"))
        {
          state_index = 0;
          show_state_or_transition_menu(simulation, state_index);
        }
        else if ((s == "u") || (s == "undo"))
        {
          handle_undo(simulation, state_index);
        }
        else if ((s == "r") || (s == "redo"))
        {
          handle_redo(simulation, state_index);
        }
        else if ((s.substr(0,2) == "g ") || (s.substr(0,5) == "goto "))
        {
          handle_goto(s, simulation, state_index);
        }
        else if ((s == "t") || (s == "trace"))
        {
          handle_trace(simulation, state_index);
        }
        else if ((s.substr(0,2) == "s ") || (s.substr(0,5) == "save "))
        {
          handle_save(s, simulation);
        }
        else if ((s.substr(0,2) == "l ") || (s.substr(0,5) == "load "))
        {
          handle_load(s, simulation, state_index);
        }
        else
        {
          std::cout << "unknown command (try 'h' for help)" << std::endl;
        }
        
      }

      return true;
    }
};

int main(int argc, char** argv)
{
  return sim_tool().execute(argc, argv);
}
