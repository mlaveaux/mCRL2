// Author(s): anonymous, Thomas Neele
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include "mcrl2/atermpp/aterm_io_binary.h"
#include "mcrl2/core/load_aterm.h"
#include "mcrl2/data/data_io.h"
#include "mcrl2/pbes/algorithms.h"
#include "mcrl2/pbes/detail/pbes_io.h"
#include "mcrl2/pbes/io.h"
#include "mcrl2/pbes/extended_pbes.h"
#include "mcrl2/pbes/join.h"
#include "mcrl2/pbes/parse.h"
#include "mcrl2/pbes/normal_forms.h"



namespace mcrl2::pbes_system
{

const std::vector<utilities::file_format>& pbes_file_formats()
{
  static std::vector<utilities::file_format> result;
  if (result.empty())
  {
    result.emplace_back("pbes", "PBES in internal format", false);
    result.back().add_extension("pbes");
    result.emplace_back("text", "PBES in textual (mCRL2) format", true);
    result.back().add_extension("txt");
    result.emplace_back("bes", "BES in internal format", false);
    result.back().add_extension("bes");
    result.emplace_back("pgsolver", "BES in PGSolver format", true);
    result.back().add_extension("gm");
    result.back().add_extension("pg");
  }
  return result;
}

void save_pbes(const pbes& pbes,
               std::ostream& stream,
               utilities::file_format format)
{
  if (format == utilities::file_format())
  {
    format = pbes_format_internal();
  }
  mCRL2log(log::verbose) << "Saving result in " << format.shortname() << " format..." << std::endl;
  if (format == pbes_format_internal() || (format == pbes_format_internal_bes() && pbes_system::algorithms::is_bes(pbes)))
  {
    atermpp::binary_aterm_ostream(stream) << pbes;
  }
  else if (format == pbes_format_pgsolver() && pbes_system::algorithms::is_bes(pbes))
  {
    save_bes_pgsolver(pbes, stream);
  }
  else if (format == pbes_format_text())
  {
    stream << pbes;
  }
  else
  {
    throw mcrl2::runtime_error("Trying to save PBES in non-PBES format (" + format.shortname() + ")");
  }
}

void load_pbes(pbes& pbes, std::istream& stream, utilities::file_format format, const std::string& /*source*/)
{
  if (format == utilities::file_format())
  {
    format = pbes_format_internal();
  }
  mCRL2log(log::verbose) << "Loading PBES in " << format.shortname() << " format..." << std::endl;
  if (format == pbes_format_internal() || format == pbes_format_internal_bes())
  {
    atermpp::binary_aterm_istream(stream) >> pbes;
  }
  else
  if (format == pbes_format_text())
  {
    throw mcrl2::runtime_error("Reading extended PBES from text format is not supported.");
  }
  else
  {
    throw mcrl2::runtime_error("Trying to load PBES from non-PBES format (" + format.shortname() + ")");
  }
}

void save_pbes(const pbes& pbes, const std::string& filename,
               utilities::file_format format,
               bool welltypedness_check)
{
  if (welltypedness_check)
  {
    assert(pbes.is_well_typed());
  }
  if (format == utilities::file_format())
  {
    format = guess_format(filename);
  }

  if (filename.empty() || filename == "-")
  {
    save_pbes(pbes, std::cout, format);
  }
  else
  {
    std::ofstream filestream(filename,(format.text_format()?std::ios_base::out: std::ios_base::binary));
    if (!filestream.good())
    {
      throw mcrl2::runtime_error("Could not open file " + filename);
    }
    save_pbes(pbes, filestream, format);
  }
}

void load_pbes(pbes& pbes,
               const std::string& filename,
               utilities::file_format format)
{
  if (format == utilities::file_format())
  {
    format = guess_format(filename);
  }
  if (filename.empty() || filename == "-")
  {
    load_pbes(pbes, std::cin, format);
  }
  else
  {
    std::ifstream filestream(filename,(format.text_format()?std::ios_base::in: std::ios_base::binary));
    if (!filestream.good())
    {
      throw mcrl2::runtime_error("Could not open file " + filename);
    }
    load_pbes(pbes, filestream, format, core::detail::file_source(filename));
  }
}

void save_extended_pbes(const extended_pbes& pbes,
               std::ostream& stream,
               utilities::file_format format)
{
  if (format == utilities::file_format())
  {
    format = pbes_format_internal();
  }
  mCRL2log(log::verbose) << "Saving result in " << format.shortname() << " format..." << std::endl;
  if (format == pbes_format_internal() || (format == pbes_format_internal_bes() && pbes_system::algorithms::is_bes(pbes.transformed_pbes)))
  {
    atermpp::binary_aterm_ostream(stream) << pbes;
  }
  else if (format == pbes_format_pgsolver() && pbes_system::algorithms::is_bes(pbes.transformed_pbes))
  {
    save_bes_pgsolver(pbes.transformed_pbes, stream);
  }
  else if (format == pbes_format_text())
  {
    stream << print(pbes, true);
  }
  else
  {
    throw mcrl2::runtime_error("Trying to save PBES in non-PBES format (" + format.shortname() + ")");
  }
}

void load_extended_pbes(extended_pbes& pbes,
               const std::string& filename,
               utilities::file_format format)
{
  if (format == utilities::file_format())
  {
    format = guess_format(filename);
  }
  if (filename.empty() || filename == "-")
  {
    load_extended_pbes(pbes, std::cin, format);
  }
  else
  {
    std::ifstream filestream(filename,(format.text_format()?std::ios_base::in: std::ios_base::binary));
    if (!filestream.good())
    {
      throw mcrl2::runtime_error("Could not open file " + filename);
    }
    load_extended_pbes(pbes, filestream, format, core::detail::file_source(filename));
  }
}

void load_extended_pbes(extended_pbes& pbes, std::istream& stream, utilities::file_format format, const std::string& /*source*/)
{
  if (format == utilities::file_format())
  {
    format = pbes_format_internal();
  }
  mCRL2log(log::verbose) << "Loading PBES in " << format.shortname() << " format..." << std::endl;
  if (format == pbes_format_internal() || format == pbes_format_internal_bes())
  {
    atermpp::binary_aterm_istream(stream) >> pbes;
  }
  else
  if (format == pbes_format_text())
  {
    throw mcrl2::runtime_error("Reading extended PBES from text format is not supported.");
  }
  else
  {
    throw mcrl2::runtime_error("Trying to load PBES from non-PBES format (" + format.shortname() + ")");
  }
}

void save_extended_pbes(const extended_pbes& pbes, const std::string& filename,
               utilities::file_format format,
               bool welltypedness_check)
{
  if (welltypedness_check)
  {
    // assert(pbes.is_well_typed()); FIXME
  }
  if (format == utilities::file_format())
  {
    format = guess_format(filename);
  }

  if (filename.empty() || filename == "-")
  {
    save_extended_pbes(pbes, std::cout, format);
  }
  else
  {
    std::ofstream filestream(filename,(format.text_format()?std::ios_base::out: std::ios_base::binary));
    if (!filestream.good())
    {
      throw mcrl2::runtime_error("Could not open file " + filename);
    }
    save_extended_pbes(pbes, filestream, format);
  }
}

// transforms OpId to OpIdNoIndex
static atermpp::aterm remove_index_impl(const atermpp::aterm& x)
{
  if (x.function() == core::detail::function_symbol_OpId())
  {
    return atermpp::aterm(core::detail::function_symbol_OpIdNoIndex(), x.begin(), --x.end());
  }
  return x;
}

// transforms OpIdNoIndex to OpId
static atermpp::aterm add_index_impl(const atermpp::aterm& x)
{
  if (x.function() == core::detail::function_symbol_OpIdNoIndex())
  {
    const data::function_symbol& y = reinterpret_cast<const data::function_symbol&>(x);
    return data::function_symbol(y.name(), y.sort());
  }
  return x;
}

inline atermpp::aterm_ostream& operator<<(atermpp::aterm_ostream& stream, const pbes_equation& equation)
{
  stream << equation.symbol();
  stream << equation.variable();
  stream << equation.formula();
  return stream;
}

inline atermpp::aterm_istream& operator>>(atermpp::aterm_istream& stream, pbes_equation& equation)
{
  fixpoint_symbol symbol;
  propositional_variable var;
  pbes_expression expression;

  stream >> symbol;
  stream >> var;
  stream >> expression;

  equation = pbes_equation(symbol, var, expression);

  return stream;
}

atermpp::aterm pbes_marker()
{
  return atermpp::aterm(atermpp::function_symbol("parameterised_boolean_equation_system", 0));
}

atermpp::aterm_ostream& operator<<(atermpp::aterm_ostream& stream, const pbes& pbes)
{
  atermpp::aterm_stream_state state(stream);
  stream << remove_index_impl;

  stream << pbes_marker();
  stream << pbes.data();
  stream << pbes.global_variables();
  stream << pbes.equations();
  stream << pbes.initial_state();
  return stream;
}

atermpp::aterm_istream& operator>>(atermpp::aterm_istream& stream, pbes& pbes)
{
  atermpp::aterm_stream_state state(stream);
  stream >> add_index_impl;

  try
  {
    atermpp::aterm marker;
    stream >> marker;

    if (marker != pbes_marker())
    {
      throw mcrl2::runtime_error("Stream does not contain a parameterised boolean equation system (PBES).");
    }

    data::data_specification data;
    std::set<data::variable> global_variables;
    std::vector<pbes_equation> equations;
    propositional_variable_instantiation initial_state;

    stream >> data;
    stream >> global_variables;
    stream >> equations;
    stream >> initial_state;

    pbes = pbes_system::pbes(data, global_variables, equations, initial_state);

    // Add all the sorts that are used in the specification
    // to the data specification. This is important for those
    // sorts that are built in, because these are not explicitly
    // declared.
    complete_data_specification(pbes);
  }
  catch (std::exception& ex)
  {
    mCRL2log(log::error) << ex.what() << "\n";
    throw mcrl2::runtime_error(std::string("Error reading parameterised boolean equation system (PBES)."));
  }

  return stream;
}

namespace detail
{

pbes load_pbes(const std::string& filename)
{
  pbes result;
  if (filename.empty() || filename == "-")
  {
    atermpp::binary_aterm_istream(std::cin) >> result;
  }
  else
  {
    std::ifstream from(filename, std::ifstream::in | std::ifstream::binary);
    atermpp::binary_aterm_istream(from) >> result;
  }
  return result;
}

void save_pbes(const pbes& pbesspec, const std::string& filename)
{
  if (filename.empty() || filename == "-")
  {
    atermpp::binary_aterm_ostream(std::cout) << pbesspec;
  }
  else
  {
    std::ofstream to(filename, std::ofstream::out | std::ofstream::binary);
    if (!to.good())
    {
      throw mcrl2::runtime_error("Could not write to filename " + filename);
    }
    atermpp::binary_aterm_ostream(to) << pbesspec;
  }
}

extended_pbes load_extended_pbes(const std::string& filename)
{
  extended_pbes result;
  if (filename.empty() || filename == "-")
  {
    atermpp::binary_aterm_istream(std::cin) >> result;
  }
  else
  {
    std::ifstream from(filename, std::ifstream::in | std::ifstream::binary);
    atermpp::binary_aterm_istream(from) >> result;
  }
  return result;
}

void save_extended_pbes(const extended_pbes& pbesspec, const std::string& filename)
{
  if (filename.empty() || filename == "-")
  {
    atermpp::binary_aterm_ostream(std::cout) << pbesspec;
  }
  else
  {
    std::ofstream to(filename, std::ofstream::out | std::ofstream::binary);
    if (!to.good())
    {
      throw mcrl2::runtime_error("Could not write to filename " + filename);
    }
    atermpp::binary_aterm_ostream(to) << pbesspec;
  }
}

} // namespace detail

/// \brief Conversion to atermappl.
/// \return The PBES converted to aterm format.
atermpp::aterm pbes_to_aterm(const pbes& p)
{
  atermpp::aterm global_variables = atermpp::aterm(core::detail::function_symbol_GlobVarSpec(),
                                                             data::variable_list(p.global_variables().begin(),
                                                                                 p.global_variables().end()));

  atermpp::aterm_list eqn_list;
  const std::vector<pbes_equation>& eqn = p.equations();
  for (auto i = eqn.rbegin(); i != eqn.rend(); ++i)
  {
    atermpp::aterm a = pbes_equation_to_aterm(*i);
    eqn_list.push_front(a);
  }
  atermpp::aterm equations = atermpp::aterm(core::detail::function_symbol_PBEqnSpec(), eqn_list);
  atermpp::aterm initial_state = atermpp::aterm(core::detail::function_symbol_PBInit(), p.initial_state());
  atermpp::aterm result;

  result = atermpp::aterm(core::detail::function_symbol_PBES(),
             data::detail::data_specification_to_aterm(p.data()),
             global_variables,
             equations,
             initial_state);

  return result;
}

inline
atermpp::aterm extended_pbes_marker()
{
  return atermpp::aterm(atermpp::function_symbol("extended_parameterised_boolean_equation_system", 0));
}

inline
atermpp::aterm_ostream& operator<<(atermpp::aterm_ostream& stream, const extended_pbes& pbes)
{
  stream << extended_pbes_marker();
  stream << pbes.transformed_pbes;
  stream << pbes.source_pbes;
  stream << pbes.source_lps;
  stream << pbes.transformations;
  return stream;
}

inline
atermpp::aterm_istream& operator>>(atermpp::aterm_istream& stream, extended_pbes& pbes)
{  
  atermpp::aterm marker;
  stream >> marker;

  if (marker == pbes_marker())
  {
    // Deal with legacy files
    stream >> pbes.transformed_pbes;
    return stream;
  }

  if (marker != extended_pbes_marker())
  {
    throw mcrl2::runtime_error("Stream does not contain an extended parameterised boolean equation system (PBES).");
  }
  
  stream >> pbes.transformed_pbes;
  stream >> pbes.source_pbes;
  stream >> pbes.source_lps;
  stream >> pbes.transformations;
  return stream;
}

} // namespace mcrl2::pbes_system


