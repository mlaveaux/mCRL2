// Author(s): Muck van Weerdenburg, Jan Friso Groote
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

/** \file
 *
 * \brief This include file contains routines to read and write
 *        labelled transitions from and to files and from streams.
 * \details This is the LTS library's main header file. It declares all publicly
 * accessible data structures that form the interface of the library.
 * \author Muck van Weerdenburg
 */

#ifndef MCRL2_LTS_LTS_IO_H
#define MCRL2_LTS_LTS_IO_H

#include <iostream>

#include "mcrl2/lps/io.h"
#include "mcrl2/lts/detail/lts_convert.h"

#include "mcrl2/atermpp/aterm_io.h"

namespace mcrl2::lts
{

/// \brief Type for data files that contain extra information for an lts in .aut or .fsm
///         format. Typically this is a data_specification (data_e), a linear process (lps_e) or
///         an .mcrl2 file. The value none_e indicates that no information is available.
using data_file_type_t = enum
{
  none_e,
  data_e,
  lps_e,
  mcrl2_e
};

namespace detail
{

/** \brief Determines the LTS format from a filename by its extension.
 * \param[in] s The name of the file of which the format will be determined.
 * \param[in] be_verbose If true, messages about the detected file format are printed in verbose mode. 
 * \return The LTS format based on the extension of \a s.
 * If no supported format can be determined from the extension then \a
 * lts_none is returned.  */
lts_type guess_format(std::string const& s, const bool be_verbose=true);

/** \brief Determines the LTS format from a format specification string.
 * \details This can be any of the following strings:
 * \li "aut" for the Ald&eacute;baran format;
 * \li "fsm" for the FSM format;
 * \li "dot" for the GraphViz format;
 *
 * \param[in] s The format specification string.
 * \return The LTS format based on the value of \a s.
 * If no supported format can be determined then \a lts_none is returned.
 */
lts_type parse_format(std::string const& s);

/** \brief Gives a string representation of an LTS format.
 * \details This is the "inverse" of \ref parse_format.
 * \param[in] type The LTS format.
 * \return The name of the LTS format specified by \a type. */
std::string string_for_type(const lts_type type);

/** \brief Gives the filename extension associated with an LTS format.
 * \param[in] type The LTS format.
 * \return The filename extension of the LTS format specified by \a type.
 */
std::string extension_for_type(const lts_type type);

/** \brief Gives the MIME type associated with an LTS format.
 * \param[in] type The LTS format.
 * \return The MIME type of the LTS format specified by \a type.
 */
std::string mime_type_for_type(const lts_type type);

/** \brief Gives the set of all supported LTS formats.
 * \return The set of all supported LTS formats. */
const std::set<lts_type>& supported_lts_formats();

/** \brief Gives a textual list describing supported LTS formats.
 * \param[in] default_format The format that should be marked as default
 *                           (or \a lts_none for no default).
 * \param[in] supported      The formats that should be considered
 *                           supported.
 * \return                   A string containing lines of the form
 *                           "  'name' for the ... format". Every line
 *                           except the last is terminated with '\\n'. */
std::string supported_lts_formats_text(lts_type default_format = lts_none, const std::set<lts_type>& supported = supported_lts_formats());

/** \brief Gives a textual list describing supported LTS formats.
 * \param[in] supported      The formats that should be considered
 *                           supported.
 * \return                   A string containing lines of the form
 *                           "  'name' for the ... format". Every line
 *                           except the last is terminated with '\\n'. */
std::string supported_lts_formats_text(const std::set<lts_type>& supported);

/** \brief Gives a list of extensions for supported LTS formats.
 * \param[in] sep       The separator to use between each extension.
 * \param[in] supported The formats that should be considered supported.
 * \return              A string containing a list of extensions of the
 *                      formats in \a supported, separated by \a sep.
 *                      E.g. "*.aut,*.lts" */
std::string lts_extensions_as_string(const std::string& sep = ",", const std::set<lts_type>& supported = supported_lts_formats());

/** \brief Gives a list of extensions for supported LTS formats.
 * \param[in] supported The formats that should be considered supported.
 * \return              A string containing a list of extensions of the
 *                      formats in \a supported, separated by \a ','.
 *                      E.g. "*.aut,*.lts" */
std::string lts_extensions_as_string(const std::set<lts_type>& supported);


inline void read_lps_context(const std::string& data_file,
                             data::data_specification& data,
                             process::action_label_list& action_labels,
                             data::variable_list& process_parameters)
{
  lps::stochastic_specification spec;
  load_lps(spec, data_file);
  data = spec.data();
  action_labels = spec.action_labels();
  process_parameters = spec.process().process_parameters();
}

inline void read_data_context(const std::string& data_file,
                              data::data_specification& data,
                              process::action_label_list& action_labels)
{
  // Add "init delta;" to the file to make the input a proper data specification.
  // "init delta;" must appear at the front, such that a second init clause will be
  // reported as wrong. 
  const std::string input="init delta; " + utilities::read_text(data_file);
  process::process_specification procspec =
      process::parse_process_specification(input);
  data=procspec.data();
  action_labels=procspec.action_labels();
}

inline void read_mcrl2_context(const std::string& data_file,
                               data::data_specification& data,
                               process::action_label_list& action_labels)
{
  process::process_specification procspec =
      process::parse_process_specification(utilities::read_text(data_file));
  data = procspec.data();
  action_labels = procspec.action_labels();
}

// converts an arbitrary lts to lts_lts_t
template <class LTS_TYPE_IN, class LTS_TYPE_OUT>
inline void convert_to_lts_lts(LTS_TYPE_IN& src, 
                               LTS_TYPE_OUT& dest, 
                               const data_file_type_t extra_data_file_type,
                               const std::string& extra_data_file_name)
{
  std::string data_file;
  data::data_specification data;
  process::action_label_list action_labels;
  data::variable_list process_parameters;
  bool extra_data_is_defined=true;
  switch (extra_data_file_type)
  {
    case data_e:  
    { 
      read_data_context(extra_data_file_name,data,action_labels); 
      break;
    }
    case lps_e:   
    {
      read_lps_context(extra_data_file_name,data,action_labels,process_parameters); 
      break;
    }
    case mcrl2_e: 
    {
      read_mcrl2_context(extra_data_file_name,data,action_labels); 
      break;
    }
    default:
    {
      extra_data_is_defined = false;
      mCRL2log(log::info) << "No data and action label specification is provided. Only the standard data types and no action labels can be used." << std::endl; break;
    }
  }
  lts_convert(src, dest, data, action_labels, process_parameters, extra_data_is_defined);
}
} //  namespace detail

/** \brief Loads an lts of the indicated type, transforms it to an lts of the form lts_lts_t using the additional data parameters.
 *  \details The file can refer to any file in lts, aut, fsm, or dot
 *           format. After reading it is is translated into .lts format. For this 
 *           a file is read with the name extra_data_file, which is interpreted 
 *           as a data specification if extra_data_file_type has type data_e, a linear process specification
 *           if it has value lps_e, and an mcrl2 file if it has value mcrl2_e.
 *  \param[out] result The lts in which the transition system is put. 
 *  \param[in] infilename The name of the file containing the lts.
 *  \param[in] type The type of the lts file, i.e. .lts, .fsm, .dot or .aut.
 *  \param[in] extra_data_file_type The type of the file containing extra information, such as a data specification.
 *  \param[in] extra_data_file_name The name of the file containing extra information. */
template <class LTS_TYPE>
inline void load_lts(LTS_TYPE& result, 
                     const std::string& infilename,
                     lts_type type,  
                     const data_file_type_t extra_data_file_type=none_e,
                     const std::string& extra_data_file_name="")
{
  switch (type)
  {
    case lts_lts:
    case lts_lts_probabilistic:
    {
      if (extra_data_file_type != none_e)
      {
        mCRL2log(log::warning) << "The lts file comes with a data specification. Ignoring the extra data and action label specification provided." << std::endl;
      }
      result.load(infilename);
      break;
    }
    case lts_none:
      mCRL2log(log::warning) << "Cannot determine type of input. Assuming .aut.\n";
      [[fallthrough]]; // For the default (lts_none) load as aut file.
    case lts_aut:
    case lts_aut_probabilistic:
    {
      probabilistic_lts_aut_t l;
      l.load(infilename);
      convert_to_lts_lts(l, result,extra_data_file_type,extra_data_file_name);
      break;
    }
    case lts_fsm:
    case lts_fsm_probabilistic:
    {
      probabilistic_lts_fsm_t l;
      l.load(infilename);
      convert_to_lts_lts(l, result,extra_data_file_type,extra_data_file_name);
      break;
    }
    case lts_dot:
    {
      throw mcrl2::runtime_error("Reading of .dot files is not supported anymore.");
    }
  }
}


/** \brief Read a labelled transition system and return it in fsm format.
 *  \details The file can refer to any file in lts, aut, fsm, or dot
 *           format. After reading it is is attempted to translate it into
 *           fsm format.
 *  \param[in] path A string with the name of the file.
 *  \param[out] l The lts in which the transition system is put. */
inline void load_lts_as_fsm_file(const std::string& path, lts_fsm_t& l)
{
  const lts_type intype = mcrl2::lts::detail::guess_format(path);
  switch (intype)
  {
    case lts_lts:
    case lts_lts_probabilistic:
    {
      lts_lts_t l1;
      l1.load(path);
      detail::lts_convert(l1,l);
      return;
    }
    case lts_none:
      mCRL2log(log::warning) << "Cannot determine type of input. Assuming .aut.\n";
      [[fallthrough]]; // For the default (lts_none) load as aut file.
    case lts_aut:
    case lts_aut_probabilistic:
    {
      lts_aut_t l1;
      l1.load(path);
      detail::lts_convert(l1,l);
      return;
    }
    case lts_fsm:
    case lts_fsm_probabilistic:
    {
      l.load(path);
      return;
    }
    case lts_dot:
    {
      throw mcrl2::runtime_error("Reading of dot files is not supported.");
    }
  }
}

// Proper interface.

/// \brief Read a (probabilistic) LTS from the given stream.
atermpp::aterm_istream& operator>>(atermpp::aterm_istream& stream, lts_lts_t& lts);
atermpp::aterm_istream& operator>>(atermpp::aterm_istream& stream, probabilistic_lts_lts_t& lts);

/// \brief Write a (probabilistic) LTS to the given stream at once.
atermpp::aterm_ostream& operator<<(atermpp::aterm_ostream& stream, const lts_lts_t& lts);
atermpp::aterm_ostream& operator<<(atermpp::aterm_ostream& stream, const probabilistic_lts_lts_t& lts);

// Streaming an LTS to disk:
//  write_lts_header(data_spec, parameters, action_labels)
//
// In any order:
//  Write transitions (to, label, from), where 'to' and 'from' are indices and 'label' the multi_action, as necessary.
//  Write state labels (state_label_lts) in their order such that writing the i-th state label belongs to state with index i.
//  Write the initial state.

/// \brief Writes the start of an LTS stream.
void write_lts_header(atermpp::aterm_ostream& stream,
  const data::data_specification& data,
  const data::variable_list& parameters,
  const process::action_label_list& action_labels);

/// \brief Write a transition to the LTS stream.
void write_transition(atermpp::aterm_ostream& stream, std::size_t from, const lps::multi_action& label, std::size_t to);
void write_transition(atermpp::aterm_ostream& stream, std::size_t from, const lps::multi_action& label, const probabilistic_lts_lts_t::probabilistic_state_t& to);

/// \brief Write a state label to the LTS stream.
void write_state_label(atermpp::aterm_ostream& stream, const state_label_lts& label);

/// \brief Write the initial state to the LTS stream.
void write_initial_state(atermpp::aterm_ostream& stream, std::size_t index);

} // namespace mcrl2::lts

#endif // MCRL2_LTS_LTS_IO_H
