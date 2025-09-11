// Author(s): Muck van Weerdenburg, Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file liblts_aut.cpp

#include <fstream>
#include "mcrl2/utilities/unordered_map.h"
#include "mcrl2/lts/lts_aut.h"
#include "mcrl2/lts/detail/liblts_swap_to_from_probabilistic_lts.h"


using namespace mcrl2::lts;

static void read_newline(std::istream& is, const std::size_t line_no)
{
  char ch;
  is.get(ch);

  // Skip over spaces
  while (ch == ' ' && is.good())
  {
    is.get(ch);
  }

  // Windows systems typically have a carriage return before a newline.
  if (ch == '\r' && is.good())
  {
    is.get(ch);
  }

  if (ch != '\n' && !is.eof()) // Last line does not need to be terminated with an eoln.
  {
    if (line_no==1)
    {
      throw mcrl2::runtime_error("Expect a newline after the header des(...,...,...).");
    }
    else
    {
      throw mcrl2::runtime_error("Expect a newline after the transition at line " + std::to_string(line_no) + ".");
    }
  }
}

// reads a number, puts it in s, and reads one extra character, which must be either a space or a closing bracket.
static void read_natural_number_to_string(std::istream& is, std::string& s, const std::size_t line_no)
{
  assert(s.empty());
  char ch;
  is >> std::skipws >> ch;
  for( ; isdigit(ch) ; is.get(ch))
  {
    s.push_back(ch);
  }
  is.putback(ch);
  if (s.empty())
  {
    throw mcrl2::runtime_error("Expect a number at line " + std::to_string(line_no) + ".");
  }
}

template <class AUT_LTS_TYPE>
static std::size_t find_label_index(const std::string& s, mcrl2::utilities::unordered_map < action_label_string, std::size_t >& labs, AUT_LTS_TYPE& l)
{
  std::size_t label;

  assert(labs.at(action_label_string::tau_action())==0);
  action_label_string as(s);
  const mcrl2::utilities::unordered_map < action_label_string, std::size_t >::const_iterator i=labs.find(as);
  if (i==labs.end())
  {
    label=l.add_action(as);
    labs[as]=label;
  }
  else
  {
    label=i->second;
  }
  return label;
}

static void check_state(std::size_t state, std::size_t number_of_states, std::size_t line_no)
{
  if (state>=number_of_states)
  {
    throw mcrl2::runtime_error("The state number " + std::to_string(state) + " is not below the number of states (" +
                               std::to_string(number_of_states) + ").  Found at line " + std::to_string(line_no) + ".");
  }
} 

static void check_states(mcrl2::lts::probabilistic_lts_aut_t::probabilistic_state_t& probability_state,
                         std::size_t number_of_states, std::size_t line_no)
{
  if (probability_state.size()<=1) // This is a simple probabilistic state. 
  {
    check_state(probability_state.get(), number_of_states, line_no);
  }
  else // The state consists of a vector of states and probability pairs. 
  {
    for(mcrl2::lts::probabilistic_lts_aut_t::probabilistic_state_t::state_probability_pair& p: probability_state)
    {
      check_state(p.state(), number_of_states, line_no);
    }
  }
} 

// This procedure tries to read states, indicated by numbers
// with in between fractions of the shape number/number. The
// last state number is put in state. The remainder as pairs
// in the vector. Typical expected input is 3 2/3 4 1/6 78 1/6 3.
static void read_probabilistic_state(
  std::istream& is,
  mcrl2::lts::probabilistic_lts_aut_t::probabilistic_state_t& result,
  const std::size_t line_no)
{
  assert(result.size()==0);

  std::size_t state;

  is >> std::skipws >> state;

  if (!is.good())
  {
    throw mcrl2::runtime_error("Expect a state number at line " + std::to_string(line_no) + ".");
  }

  // Check whether the next character is a digit. If so a probability follows.
  char ch;
  is >> std::skipws >> ch;
  is.putback(ch);

  if (!isdigit(ch))
  {
    // There is only a single state.
    result.set(state);
    return;
  }
  bool ready=false;

  mcrl2::utilities::probabilistic_arbitrary_precision_fraction remainder=mcrl2::utilities::probabilistic_arbitrary_precision_fraction::one();
  while (is.good() && !ready)
  {
    // Now read a probabilities followed by the next state.
    std::string enumerator;
    read_natural_number_to_string(is,enumerator,line_no);
    char ch;
    is >> std::skipws >> ch;
    if (ch != '/')
    {
      throw mcrl2::runtime_error("Expect a / in a probability at line " + std::to_string(line_no) + ".");
    }

    std::string denominator;
    read_natural_number_to_string(is,denominator,line_no);
    mcrl2::utilities::probabilistic_arbitrary_precision_fraction frac(enumerator,denominator);
    remainder=remainder-frac;
    result.add(state, frac);
    
    is >> std::skipws >> state;

    if (!is.good())
    {
      throw mcrl2::runtime_error("Expect a state number at line " + std::to_string(line_no) + ".");
    }

    // Check whether the next character is a digit.
    
    is >> std::skipws >> ch;
    is.putback(ch);

    if (!isdigit(ch))
    {
      ready=true;
    }
  }
  
  result.add(state, remainder);
}


static void read_aut_header(
  std::istream& is,
  mcrl2::lts::probabilistic_lts_aut_t::probabilistic_state_t& initial_state,
  std::size_t& num_transitions,
  std::size_t& num_states)
{
  std::string s;
  is.width(3);
  is >> std::skipws >> s;

  if (s!="des")
  {
    throw mcrl2::runtime_error("Expect an .aut file to start with 'des'.");
  }

  char ch;
  is >> std::skipws >> ch;

  if (ch != '(')
  {
    throw mcrl2::runtime_error("Expect an opening bracket '(' after 'des' in the first line of a .aut file.");
  }

  read_probabilistic_state(is,initial_state,1);

  is >> std::skipws >> ch;
  if (ch != ',')
  {
    throw mcrl2::runtime_error("Expect a comma after the first number in the first line of a .aut file.");
  }

  is >> std::skipws >> num_transitions;

  is >> std::skipws >> ch;
  if (ch != ',')
  {
    throw mcrl2::runtime_error("Expect a comma after the second number in the first line of a .aut file.");
  }

  is >> std::skipws >> num_states;

  is >> ch;

  if (ch != ')')
  {
    throw mcrl2::runtime_error("Expect a closing bracket ')' after the third number in the first line of a .aut file.");
  }

  read_newline(is,1);
}

static bool read_initial_part_of_an_aut_transition(
  std::istream& is,
  std::size_t& from,
  std::string& label,
  const std::size_t line_no)
{
  char ch;
  is >> std::skipws >> ch;
  if (is.eof())
  {
    return false;
  }
  if (ch == 0x04) // found EOT character that separates two files
  {
    return false;
  }

  is >> std::skipws >> from;

  is >> std::skipws >> ch;
  if (ch != ',')
  {
    throw mcrl2::runtime_error("Expect that the first number is followed by a comma at line " + std::to_string(line_no) + ".");
  }

  is >> std::skipws >> ch;
  if (ch == '"')
  {
    label="";
    // In case the label is using quotes whitespaces
    // in the label are preserved. 
    is >> std::noskipws >> ch;
    while ((ch != '"') && !is.eof())
    {
      label.push_back(ch);
      is >> ch;
    }
    if (ch != '"')
    {
      throw mcrl2::runtime_error("Expect that the second item is a quoted label (using \") at line " + std::to_string(line_no) + ".");
    }
    is >> std::skipws >> ch;
  }
  else
  {
    // In case the label is not within quotes,
    // whitespaces are removed from the label. 
    label = ch;
    is >> ch;
    while ((ch != ',') && !is.eof())
    {
      label.push_back(ch);
      is >> ch;
    }
  }

  if (ch != ',')
  {
    throw mcrl2::runtime_error("Expect a comma after the quoted label at line " + std::to_string(line_no) + ".");
  }
  
  return true;
}

static bool read_aut_transition(
  std::istream& is,
  std::size_t& from,
  std::string& label,
  mcrl2::lts::probabilistic_lts_aut_t::probabilistic_state_t& target_probabilistic_state,
  const std::size_t line_no)
{
  if (!read_initial_part_of_an_aut_transition(is,from,label,line_no))
  {
    return false;
  }

  read_probabilistic_state(is,target_probabilistic_state,line_no);

  char ch;
  is >> ch;
  if (ch != ')')
  {
    throw mcrl2::runtime_error("Expect a closing bracket at the end of the transition at line " + std::to_string(line_no) + ".");
  }

  read_newline(is,line_no);
  return true;
}

static bool read_aut_transition(
  std::istream& is,
  std::size_t& from,
  std::string& label,
  std::size_t& to,
  const std::size_t line_no)
{
  if (!read_initial_part_of_an_aut_transition(is,from,label,line_no))
  {
    return false;
  }

  is >> std::skipws >> to;

  char ch;
  is >> ch;
  if (ch != ')')
  {
    throw mcrl2::runtime_error("Expect a closing bracket at the end of the transition at line " + std::to_string(line_no) + ".");
  }

  read_newline(is,line_no);
  return true;
}

static size_t add_probablistic_state(
                    mcrl2::lts::probabilistic_lts_aut_t::probabilistic_state_t& probabilistic_state,
                    probabilistic_lts_aut_t& l,
                    mcrl2::utilities::unordered_map < std::size_t, std::size_t>& indices_of_single_probabilistic_states,
                    mcrl2::utilities::unordered_map < mcrl2::lts::probabilistic_lts_aut_t::probabilistic_state_t, std::size_t>& 
                                              indices_of_multiple_probabilistic_states) 
{
  std::size_t fresh_index = indices_of_single_probabilistic_states.size()+indices_of_multiple_probabilistic_states.size();
  std::size_t index;
  // Check whether probabilistic states exists. 
  if (probabilistic_state.size()<=1)
  {
    index = indices_of_single_probabilistic_states.insert(
                     std::pair< std::size_t, std::size_t>
                     (probabilistic_state.get(),fresh_index)).first->second;
  }
  else
  {
    assert(probabilistic_state.size()>1);
    index = indices_of_multiple_probabilistic_states.insert(
                     std::pair< mcrl2::lts::probabilistic_lts_aut_t::probabilistic_state_t, std::size_t>
                     (probabilistic_state,fresh_index)).first->second;
  }
 
  if (index==fresh_index)
  {
    std::size_t probabilistic_state_index=l.add_and_reset_probabilistic_state(probabilistic_state);
    assert(probabilistic_state_index==index);
    (void)probabilistic_state_index; // Avoid unused variable warning.
  }
  return index;
}


static void read_from_aut(probabilistic_lts_aut_t& l, std::istream& is)
{
  std::size_t line_no = 1;
  std::size_t ntrans = 0;
  std::size_t nstate = 0;

  mcrl2::lts::probabilistic_lts_aut_t::probabilistic_state_t initial_probabilistic_state;
  read_aut_header(is,initial_probabilistic_state,ntrans,nstate);

  // The two unordered maps below are used to determine a unique index for each probabilistic state.
  // Because most states consist of one probabilistic state, the unordered maps are duplicated into
  // indices_of_single_probabilistic_states and indices_of_multiple_probabilistic_states.
  // The map indices_of_single_probabilistic_states requires far less memory.
  mcrl2::utilities::unordered_map < std::size_t, std::size_t> indices_of_single_probabilistic_states;
  mcrl2::utilities::unordered_map < mcrl2::lts::probabilistic_lts_aut_t::probabilistic_state_t, std::size_t> indices_of_multiple_probabilistic_states;
  
  check_states(initial_probabilistic_state, nstate, line_no);

  if (nstate==0)
  {
    throw mcrl2::runtime_error("cannot parse AUT input that has no states; at least an initial state is required.");
  }

  l.set_num_states(nstate,false);
  l.clear_transitions(ntrans); // Reserve enough space for the transitions.
  
  mcrl2::utilities::unordered_map < action_label_string, std::size_t > action_labels;
  action_labels[action_label_string::tau_action()]=0; // A tau action is always stored at position 0.
  l.set_initial_probabilistic_state(initial_probabilistic_state); 

  mcrl2::lts::probabilistic_lts_aut_t::probabilistic_state_t probabilistic_target_state;
  std::size_t from;
  std::string s;

  while (!is.eof())
  {
    probabilistic_target_state.clear();

    line_no++;

    if (!read_aut_transition(is,from,s,probabilistic_target_state,line_no))
    {
      break; // encountered EOF or something that is not a transition
    }

    check_state(from, nstate, line_no);
    check_states(probabilistic_target_state, nstate, line_no);
    std::size_t index = add_probablistic_state(probabilistic_target_state, l, indices_of_single_probabilistic_states, indices_of_multiple_probabilistic_states);

    l.add_transition(transition(from,find_label_index(s,action_labels,l),index));
  }

  if (ntrans != l.num_transitions())
  {
    throw mcrl2::runtime_error("number of transitions read (" + std::to_string(l.num_transitions()) +
                               ") does not correspond to the number of transition given in the header (" + std::to_string(ntrans) + ").");
  }
}

static void read_from_aut(lts_aut_t& l, std::istream& is)
{
  std::size_t line_no = 1;
  std::size_t ntrans = 0;
  std::size_t nstate = 0;

  mcrl2::lts::probabilistic_lts_aut_t::probabilistic_state_t initial_probabilistic_state;
  read_aut_header(is,initial_probabilistic_state,ntrans,nstate);
  
  if (initial_probabilistic_state.size()>1)
  {
    throw mcrl2::runtime_error("Encountered an initial probability distribution while reading an non probabilistic .aut file.");
  }

  check_states(initial_probabilistic_state, nstate, line_no);

  if (nstate==0)
  {
    throw mcrl2::runtime_error("cannot parse AUT input that has no states; at least an initial state is required.");
  }

  l.set_num_states(nstate,false);
  l.clear_transitions(ntrans); // Reserve enough space for the transitions.
  
  mcrl2::utilities::unordered_map < action_label_string, std::size_t > action_labels;
  action_labels[action_label_string::tau_action()]=0; // A tau action is always stored at position 0.
  l.set_initial_state(initial_probabilistic_state.get());

  std::size_t from;
  std::size_t to;
  std::string s;
  while (!is.eof())
  {
    line_no++;

    if (!read_aut_transition(is,from,s,to,line_no))
    {
      break; // eof encountered
    }

    check_state(from, nstate, line_no);
    check_state(to, nstate, line_no);
    l.add_transition(transition(from,find_label_index(s,action_labels,l),to));
  }

  if (ntrans != l.num_transitions())
  {
    throw mcrl2::runtime_error("number of transitions read (" + std::to_string(l.num_transitions()) +
                               ") does not correspond to the number of transition given in the header (" + std::to_string(ntrans) + ").");
  }
}


static void write_probabilistic_state(const mcrl2::lts::probabilistic_lts_aut_t::probabilistic_state_t& prob_state, std::ostream& os)
{
  mcrl2::utilities::probabilistic_arbitrary_precision_fraction previous_probability;
  bool first_element=true;
  if (prob_state.size()<=1) // This is a simple probabilistic state. 
  {
    os << prob_state.get();
  }
  else // The state consists of a vector of states and probability pairs. 
  {
    for (const mcrl2::lts::probabilistic_lts_aut_t::probabilistic_state_t::state_probability_pair& p: prob_state)
    {
      if (first_element)
      {
        os << p.state();
        previous_probability=p.probability();
        first_element=false;
      }
      else
      {
        os << " " << pp(previous_probability) << " " << p.state();
        previous_probability=p.probability();
      }
    }
  }
}

static void write_to_aut(const probabilistic_lts_aut_t& l, std::ostream& os)
{
  // Do not use "endl" below to avoid flushing. Use "\n" instead.
  os << "des (";
  write_probabilistic_state(l.initial_probabilistic_state(),os);

  os << "," << l.num_transitions() << "," << l.num_states() << ")" << "\n";

  for (const transition& t: l.get_transitions())
  {
    os << "(" << t.from() << ",\"" << pp(l.action_label(l.apply_hidden_label_map(t.label()))) << "\",";
    write_probabilistic_state(l.probabilistic_state(t.to()),os);
    os << ")" << "\n";
  }
}

static void write_to_aut(const lts_aut_t& l, std::ostream& os)
{
  // Do not use "endl" below to avoid flushing. Use "\n" instead.
  os << "des (" << l.initial_state() << "," << l.num_transitions() << "," << l.num_states() << ")" << "\n"; 

  for (const transition& t: l.get_transitions())
  {
    os << "(" << t.from() << ",\"" 
       << pp(l.action_label(l.apply_hidden_label_map(t.label()))) << "\"," 
       << t.to() << ")" << "\n";
  }
}


namespace mcrl2::lts
{

void probabilistic_lts_aut_t::load(const std::string& filename)
{
  if (filename=="" || filename=="-")
  {
    read_from_aut(*this, std::cin);
  }
  else
  {
    std::ifstream is(filename.c_str());

    if (!is.is_open())
    {
      throw mcrl2::runtime_error("cannot open .aut file '" + filename + ".");
    }

    read_from_aut(*this,is);
    is.close();
  }
}

void probabilistic_lts_aut_t::load(std::istream& is)
{
  read_from_aut(*this,is);
}

void probabilistic_lts_aut_t::save(std::string const& filename) const
{
  if (filename=="" || filename=="-")
  {
    write_to_aut(*this, std::cout);
  }
  else
  {
    std::ofstream os(filename.c_str());

    if (!os.is_open())
    {
      throw mcrl2::runtime_error("cannot create .aut file '" + filename + ".");
      return;
    }
    write_to_aut(*this,os);
    os.close();
  }
}

void lts_aut_t::load(const std::string& filename)
{
  if (filename.empty() || filename=="-")
  {
    read_from_aut(*this, std::cin);
  }
  else
  {
    std::ifstream is(filename.c_str());

    if (!is.is_open())
    {
      throw mcrl2::runtime_error("cannot open .aut file '" + filename + ".");
    }

    read_from_aut(*this,is);
    is.close();
  }
}

void lts_aut_t::load(std::istream& is)
{
  read_from_aut(*this,is);
}

void lts_aut_t::save(std::string const& filename) const
{
  if (filename.empty() || filename=="-")
  {
    write_to_aut(*this, std::cout);
  }
  else
  {
    std::ofstream os(filename.c_str());

    if (!os.is_open())
    {
      throw mcrl2::runtime_error("cannot create .aut file '" + filename + ".");
      return;
    }
    write_to_aut(*this,os);
    os.close();
  }
}


}

