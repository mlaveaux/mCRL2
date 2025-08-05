// Author(s): Jeroen van der Wulp
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)


#include "mcrl2/utilities/platform.h"
#include "mcrl2/utilities/command_line_interface.h"

#include <locale>
#include <regex>

namespace mcrl2::utilities
{

template <typename Iter>
std::string substring(const std::string& s, Iter first, Iter last)
{
  assert(last <= s.end());
  return s.substr(first - s.begin(), last - first);
}

/// \cond INTERNAL
/**
 * Inserts newline characters while reading the input from left to right.
 * When a newline is inserted the next read newline is discarded.
 *  - if a newline is read the previous line was wrapped the first
 *
 * To force a new-line at some position add two successive newline characters.
 *
 * \param[in] input the input string
 * \param[in] indent used as fixed indentation after end of line
 * \param[in] width the maximum width of the text
 * \pre 0 < width and no word can be longer than width
 * \return string with newlines inserted such that the number of characters
 * between any two consecutive newlines is smaller than width
 **/
static std::string word_wrap(const std::string& input, const std::size_t width, const std::string& indent = "")
{
  std::ostringstream out;

  std::string            variable_indent = input.substr(0, input.find_first_not_of(" \t"));
  std::string::size_type space_left      = width;

  std::string::const_iterator i = input.begin();
  std::string::const_iterator word_start = i;

  while (i != input.end())
  {
    if (space_left - (i - word_start) < 1)   // line too long
    {

      out << std::endl << indent << variable_indent << substring(input, word_start, i);

      space_left = width - (i - word_start) - variable_indent.size();
      word_start = i;
    }
    else if (*i == '\n')
    {
      if (word_start != i)
      {
        out << substring(input, word_start, i);
      }

      ++i;

      // store new indent
      std::string::size_type start = i - input.begin();
      std::string::size_type end   = input.find_first_not_of(" \t", start);

      if (end != std::string::npos)
      {
        i              += end - start;
        variable_indent = input.substr(start, end - start);
      }
      else
      {
        variable_indent.clear();
      }

      // copy word and newline
      out << std::endl << indent << variable_indent;

      space_left = width - variable_indent.size();
      word_start = i;
    }
    else if (1 < space_left && (*i == ' ' || *i == '\t'))
    {
      out << substring(input, word_start, ++i);

      space_left -= i - word_start;
      word_start  = i;
    }
    else
    {
      ++i;
    }
  }

  out << substring(input, word_start, input.end());

  return out.str();
}

/**
 * \param[in] left_width the width of the left column
 * \param[in] right_width the width of the right column
 * \return formatted string that represents the option description
 **/
std::string interface_description::option_descriptor::textual_description(
  const std::size_t left_width, const std::size_t right_width) const
{
  std::ostringstream s;
  std::string        options;

  if (m_short != '\0')
  {
    options = "  -" + std::string(1, m_short);

    if (m_argument.get() != nullptr)
    {
      options += (m_argument->is_optional()) ?
                 "[" + m_argument->get_name() + "]" : m_argument->get_name();
    }

    options += ", ";
  }
  else
  {
    options = "      ";
  }

  options += "--" + m_long;

  if (m_argument.get() != nullptr)
  {
    options += ((m_argument->is_optional()) ?
                "[=" + m_argument->get_name() + "]" : "=" + m_argument->get_name());
  }

  if (options.size() < left_width)
  {
    s << options << std::string(left_width - options.size(), ' ');
  }
  else
  {
    s << options << std::endl << std::string(left_width, ' ');
  }

  s << word_wrap(m_description, right_width, std::string(left_width, ' ')) << std::endl;
  if (m_argument.get() != nullptr && m_argument->has_description())
  {
    std::vector< basic_argument::argument_description > arg_description(m_argument->get_description());
    for(std::vector< basic_argument::argument_description >::const_iterator i = arg_description.begin(); i != arg_description.end(); ++i)
    {
      std::string arg;
      if(i->get_short() != std::string())
      {
        arg += "'" + i->get_short() + "', ";
      }
      arg += "'" + i->get_long() + "' " + i->get_description();
      bool is_default = m_argument->get_default() == i->get_long();
      if(is_default)
      {
        arg += " (default)";
      }
      s << std::string(left_width+2, ' ')
        << word_wrap(arg, right_width, std::string(left_width+4, ' ')) << std::endl;
    }
  }

  return s.str();
}

std::string interface_description::option_descriptor::man_page_description() const
{

  std::ostringstream s;

  s << ".TP" << std::endl;

  if (m_short != '\0')
  {
    s << "\\fB-" << std::string(1, m_short) << "\\fR";

    if (m_argument.get() != nullptr)
    {
      if (m_argument->is_optional())
      {
        s << "[\\fI" << m_argument->get_name() << "\\fR]";
      }
      else
      {
        s << "\\fI" << m_argument->get_name() << "\\fR";
      }
    }

    s << ", ";
  }

  s << "\\fB--" << m_long << "\\fR";

  if (m_argument.get() != nullptr)
  {
    s << ((m_argument->is_optional()) ?
          "[=\\fI" + m_argument->get_name() + "\\fR]" :
          "=\\fI" + m_argument->get_name() + "\\fR");
  }

  s << std::endl
    << std::regex_replace(
         std::regex_replace(word_wrap(m_description, 80),
                            std::regex("\'"), 
                            std::string("\\&'")),
         std::regex("\\."), 
         "\\&.")
    << std::endl;

  if (m_argument.get() != nullptr && m_argument->has_description())
  {
    std::vector< basic_argument::argument_description > arg_description(m_argument->get_description());
    for(std::vector< basic_argument::argument_description >::const_iterator i = arg_description.begin(); i != arg_description.end(); ++i)
    {
      std::string arg;
      if(i->get_short() != std::string())
      {
        arg += "'" + i->get_short() + "', ";
      }
      arg += "'" + i->get_long() + "' " + i->get_description();
      bool is_default = m_argument->get_default() == i->get_long();
      if(is_default)
      {
        arg += " (default)";
      }

      s << std::regex_replace(
            std::regex_replace(word_wrap(arg, 80),
                               std::regex("\'"), 
                               std::string("\\&'")),
            std::regex("\\."), 
            std::string("\\&."))
       << std::endl;
    }
  }

  return s.str();
}

std::ostream& interface_description::option_descriptor::xml_page_description(std::ostream& s, const bool is_standard, unsigned int indentation) const
{
  s << std::string(indentation++, ' ') << "<option standard=\"" << (is_standard?"yes":"no") << "\">" << std::endl;

  if (m_short != '\0')
  {
    s << std::string(indentation, ' ') << "<short>" << m_short << "</short>" << std::endl;
  }

  s << std::string(indentation, ' ') << "<long>" << m_long << "</long>" << std::endl;

  s << std::string(indentation++, ' ') << "<description>";
  // Produce output line by line, such that indentation provided in the description
  // can be preserved.
  std::vector<std::string> lines = mcrl2::utilities::split(m_description, "\n");
  for (std::vector<std::string>::const_iterator i = lines.begin(); i != lines.end(); ++i)
  {
    s << *i << "<br/>" << std::endl;
  }
  s << std::string(--indentation, ' ') << "</description>" << std::endl;

  if (m_argument.get() != nullptr)
  {
    s << std::string(indentation++, ' ') << "<option_argument optional=\"" << (m_argument->is_optional()?"yes":"no") <<"\" type=\"" << (m_argument->get_type()) << "\">" << std::endl;
    s << std::string(indentation, ' ') << "<name>" << m_argument->get_name() << "</name>" << std::endl;

    if(m_argument->has_description())
    {
      s << std::string(indentation++, ' ') << "<values>" << std::endl;
      std::vector< basic_argument::argument_description > arg_description(m_argument->get_description());
      for(std::vector< basic_argument::argument_description >::const_iterator i = arg_description.begin(); i != arg_description.end(); ++i)
      {
        bool is_default = m_argument->get_default() == i->get_long();
        s << std::string(indentation++, ' ') << "<value default=\"" << (is_default?"yes":"no") << "\">" << std::endl;
        if(i->get_short() != std::string())
        {
          s << std::string(indentation, ' ') << "<short>" << i->get_short() << "</short>" << std::endl;
        }
        s << std::string(indentation, ' ') << "<long>" << i->get_long() << "</long>" << std::endl;
        s << std::string(indentation, ' ') << "<description>" << i->get_description() << "</description>" << std::endl;
        s << std::string(--indentation, ' ') << "</value>" << std::endl;
      }
      s << std::string(--indentation, ' ') << "</values>" << std::endl;
    }

    s << std::string(--indentation, ' ') << "</option_argument>" << std::endl;
  }

  s << std::string(--indentation, ' ') << "</option>" << std::endl;

  return s;
}

interface_description::option_descriptor const& interface_description::find_option(std::string const& long_identifier) const
{
  option_map::const_iterator i(m_options.find(long_identifier));

  if (i == m_options.end())
  {
    throw std::logic_error("Find operation for invalid option `" +  long_identifier + "'\n");
  }

  return i->second;
}
/// \endcond

interface_description::mandatory_argument< std::string >
make_mandatory_argument(std::string const& name, std::string const& default_value)
{

  return interface_description::mandatory_argument< std::string >(name, default_value);
}

interface_description::mandatory_argument< std::string >
make_mandatory_argument(std::string const& name)
{

  return interface_description::mandatory_argument< std::string >(name);
}

interface_description::optional_argument< std::string >
make_optional_argument(std::string const& name, std::string const& default_value)
{

  return interface_description::optional_argument< std::string >(name, default_value);
}

interface_description& interface_description::get_standard_description()
{
  static interface_description d;

  return d;
}

interface_description::interface_description(std::string const& path,
    std::string const& name, std::string const& authors,
    std::string const& what_is, std::string const& synopsis,
    std::string const& description, std::string const& known_issues) :
  m_options(get_standard_description().m_options),
  m_path(path), m_name(name), m_authors(authors),
  m_what_is(what_is), m_usage(synopsis), m_description(description),
  m_known_issues(known_issues),
  m_short_to_long(get_standard_description().m_short_to_long)
{

  m_usage = m_usage.substr(0, m_usage.find_last_of('\n'));

  // Add mandatory options
  add_hidden_option("help", "display help information", 'h');
  add_hidden_option("help-all", "display help information, including hidden and experimental options");
  add_hidden_option("version", "display version information");
  add_hidden_option("quiet", "do not display warning messages", 'q');
  add_hidden_option("verbose", "display short log messages", 'v');
  add_hidden_option("debug", "display detailed log messages", 'd');
  add_hidden_option("log-level", make_mandatory_argument<std::string>("LEVEL", ""),
                    "display log messages up to and including level; either warn, verbose, debug or trace");
}

std::string interface_description::copyright_message()
{
  return "Copyright (c) " + copyright_period() +
         " Technische Universiteit Eindhoven.\n"
         "This is free software.  You may redistribute copies of it under the\n"
         "terms of the Boost Software License <http://www.boost.org/LICENSE_1_0.txt>.\n"
         "There is NO WARRANTY, to the extent permitted by law.\n";
}

std::string interface_description::version_information() const
{
  return m_name + " mCRL2 toolset " + get_toolset_version() + "\n" +
         copyright_message() + "\nWritten by " + m_authors + ".\n";
}

std::string interface_description::textual_description(bool show_hidden) const
{
  std::ostringstream s;

  s << "Usage: " << m_path << " " << m_usage << std::endl
    << word_wrap(m_description, 80) << std::endl << std::endl;

  if (5 < m_options.size())   // tool-specific options
  {
    option_map::const_iterator        i = m_options.begin();
    short_to_long_map::const_iterator j = m_short_to_long.begin();

    s << "Options:" << std::endl;

    while (true)
    {
      option_descriptor const* option;

      // skip options without short identifier
      while (i != m_options.end() && i->second.m_short != '\0')
      {
        ++i;
      }

      if (j != m_short_to_long.end())
      {
        if (i != m_options.end())
        {
          if (i->second.m_long[0] < j->first)
          {
            option = &(i++)->second;
          }
          else
          {
            if (i->first == j->second)
            {
              ++i;
            }

            option = &m_options.find((j++)->second)->second;
          }
        }
        else
        {
          option = &m_options.find((j++)->second)->second;
        }
      }
      else if (i != m_options.end())
      {
        option = &(i++)->second;
      }
      else
      {
        break;
      }

      if ((show_hidden && option->get_long() != "quiet"
          && option->get_long() != "verbose" && option->get_long() != "debug"
          && option->get_long() != "log-level" && option->get_long() != "help"
          && option->get_long() != "help-all" && option->get_long() != "version")
          || option->m_show)
      {
        s << option->textual_description(27, 53);
      }
    }

    s << std::endl;
  }

  s << "Standard options:" << std::endl
    << m_options.find("quiet")->second.textual_description(27, 53)
    << m_options.find("verbose")->second.textual_description(27, 53)
    << m_options.find("debug")->second.textual_description(27, 53)
    << m_options.find("log-level")->second.textual_description(27, 53)
    << m_options.find("help")->second.textual_description(27, 53)
    << m_options.find("version")->second.textual_description(27, 53)
    << m_options.find("help-all")->second.textual_description(27, 53)
    << std::endl;

  if (!m_known_issues.empty())
  {
    s << "Known Issues:" << std::endl
      << word_wrap(m_known_issues, 80) << std::endl << std::endl;
  }

  s << "Report bugs at <http://www.mcrl2.org/issuetracker>." << std::endl
    << std::endl
    << "See also the manual at <http://www.mcrl2.org/web/user_manual/tools/release/" << m_name << ".html>.\n";

  return s.str();
}

inline static std::string mark_name_in_usage(std::string const& usage, const std::string& begin, const std::string& end)
{
  std::string result;
  bool        name_character = false;

  result.reserve(2 * usage.size());

  for (char i : usage)
  {
    if (i == '[')
    {
      name_character = true;
      result.append("[" + begin);
    }
    else if ((i == ' ') || (i == ']'))
    {
      if (name_character)
      {
        result.append(end);

        name_character = false;
      }
      result.append(1, i);
    }
    else
    {
      result.append(1, i);
    }
  }

  return result;
}

std::string interface_description::man_page() const
{
  std::ostringstream s;

  s.imbue(std::locale(s.getloc()));

  s << ".\\\" " << "Manual page for " << m_name << " version " << get_toolset_version() << "." << " .\\\"" << std::endl
    << ".\\\" " << "Generated from " << m_name << " --generate-man-page." << " .\\\""<< std::endl;

  // Determine month and year, to prevent using boost date/time
  time_t rawtime;
  struct tm* timeinfo;
  char buffer [80];
  time(&rawtime);
  timeinfo = localtime(&rawtime);
  strftime(buffer,80,"%B %Y",timeinfo);

  std::string name_upper = m_name;
  std::transform(name_upper.begin(), name_upper.end(), name_upper.begin(), ::toupper);
  s << ".TH " << name_upper << R"( "1" ")" << std::string(buffer) << "\" \"" << m_name << " mCRL2 toolset "
    << get_toolset_version() << R"(" "User Commands")" << std::endl;

  s << ".SH NAME" << std::endl
    << m_name << " \\- " << m_what_is << std::endl;

  s << ".SH SYNOPSIS" << std::endl
    << ".B " << m_name << std::endl
    << mark_name_in_usage(m_usage, "\\fI", "\\fR") << std::endl;

  s << ".SH DESCRIPTION" << std::endl
    << word_wrap(m_description, 80) << std::endl;

  if (0 < m_options.size())
  {
    s << ".SH OPTIONS" << std::endl
      << ".TP" << std::endl
      << "\\fIOPTION\\fR can be any of the following:" << std::endl;

    for (const auto & m_option : m_options)
    {
      option_descriptor const& option(m_option.second);

      if (option.m_show)
      {
        s << option.man_page_description();
      }
    }
  }

  s << ".TP" << std::endl;

  if (0 < m_options.size())
  {
    s << "Standard options:" << std::endl;
  }
  else
  {
    s << "\\fIOPTION\\fR can be any of the following standard options:\n";
  }
  s << m_options.find("quiet")->second.man_page_description()
    << m_options.find("verbose")->second.man_page_description()
    << m_options.find("debug")->second.man_page_description()
    << m_options.find("log-level")->second.man_page_description()
    << m_options.find("help")->second.man_page_description()
    << m_options.find("version")->second.man_page_description()
    << m_options.find("help-all")->second.man_page_description()
    << std::endl;

  if (!m_known_issues.empty())
  {
    s << ".SH \"KNOWN ISSUES\"" << std::endl
      << word_wrap(m_known_issues, 80) << std::endl;
  }

  s << ".SH AUTHOR" << std::endl
    << "Written by " << m_authors << "." << std::endl;
  s << ".SH \"REPORTING BUGS\"" << std::endl
    << "Report bugs at <http://www.mcrl2.org/issuetracker>." << std::endl;
  s << ".SH COPYRIGHT" << std::endl
    << "Copyright \\(co " + copyright_period() +
    " Technische Universiteit Eindhoven.\n"
    << ".br" << std::endl
    << "This is free software.  You may redistribute copies of it under the\n"
    "terms of the Boost Software License <http://www.boost.org/LICENSE_1_0.txt>.\n"
    "There is NO WARRANTY, to the extent permitted by law.\n";
  s << ".SH \"SEE ALSO\"" << std::endl
    << "See also the manual at <http://www.mcrl2.org/web/user_manual/tools/release/" << m_name << ".html>.\n";

  return s.str();
}

std::map<std::string, std::string> interface_description::get_long_argument_with_description()
{
  std::map<std::string, std::string> result;

  for (option_map::const_iterator i = m_options.begin(); i != m_options.end(); ++i)
  {
    option_descriptor const& option(i->second);

    result.insert(std::pair<std::string, std::string>(option.m_long ,option.m_description));
  }
  return result;
}

std::ostream& interface_description::xml_page(std::ostream& s) const
{
  unsigned int indentation = 0;

  s << std::string(indentation++, ' ') << "<tool>" << std::endl;
  s << std::string(indentation, ' ')   << "<name>" << m_name << "</name>" << std::endl;
  s << std::string(indentation, ' ')   << "<usage>" << m_usage << "</usage>" << std::endl;
  s << std::string(indentation++, ' ') << "<description>" << std::endl;

  std::vector<std::string> lines = mcrl2::utilities::split(m_description, "\n");
  for (std::vector<std::string>::const_iterator i = lines.begin(); i != lines.end(); ++i)
  {
    s << *i << "<br/>" << std::endl;
  }
  s << std::string(--indentation, ' ') << "</description>" << std::endl;

  if (0 < m_options.size())
  {
    s << std::string(indentation++, ' ') << "<options>" << std::endl;

    for (const auto & m_option : m_options)
    {
      option_descriptor const& option(m_option.second);

      if (option.m_show)
      {
        option.xml_page_description(s, false, indentation);
      }
    }
  }

  m_options.find("quiet")->second.xml_page_description(s, true, indentation);
  m_options.find("verbose")->second.xml_page_description(s, true, indentation);
  m_options.find("debug")->second.xml_page_description(s, true, indentation);
  m_options.find("log-level")->second.xml_page_description(s, true, indentation);
  m_options.find("help")->second.xml_page_description(s, true, indentation);
  m_options.find("version")->second.xml_page_description(s, true, indentation);
  m_options.find("help-all")->second.xml_page_description(s, true, indentation);

  s << std::string(--indentation, ' ') << "</options>" << std::endl;

  if (0 < m_known_issues.size())
  {
    s << std::string(indentation, ' ') << "<known_issues>" << m_known_issues << "</known_issues>" << std::endl;
  }
  s << std::string(indentation, ' ') << "<author>" << m_authors << "</author>" << std::endl;
  s << std::string(--indentation, ' ') << "</tool>" << std::endl;

  assert(indentation == 0);

  return s;
}

/**
 * The interface description specifies the available options and their
 * possible arguments. This procedure traverses the command-line arguments
 * in a left-to-right fashion and assigns those arguments to options. The
 * output is a pair of a map and a vector: the option map and the arguments
 * vector. The option map (a multimap) maps an options found as argument to
 * its argument or an empty string (in case of no argument). The arguments
 * vector contains command line arguments that were not recognised as
 * option or argument to option.
 *
 * \param[in] d interface description
 * \param[in] arguments vector (as ordered list) with command line arguments
 * \throw std::runtime_error with parse error message
 **/
void command_line_parser::collect(interface_description& d, std::vector< std::string > const& arguments)
{
  std::vector< std::string >::const_iterator i = arguments.begin();

  if (arguments.size() == 0)
  {
    return;
  }

#ifdef MCRL2_PLATFORM_MAC
  if (1 < arguments.size() && arguments[1].substr(0, 7) == "-psn_0_")   // Skip Apple serial number
  {
    ++i;
  }
#endif

  while (++i != arguments.end())
  {
    std::string const& argument(*i);

    if (1 < argument.size() && argument[0] == '-')
    {
      if (argument[1] == '-')
      {
        // Assume that the argument is a long option
        std::string option(argument, 2);

        if (option.find_first_of('=') < option.size())
        {
          // remove argument to single out the long option
          option.resize(option.find_first_of('='));
        }

        if (d.m_options.find(option) == d.m_options.end())
        {
          if (argument == "--generate-man-page" || argument == "--generate-xml")
          {
            // Special option
            m_options.insert(std::make_pair(argument.substr(2), ""));
          }
          else
          {
            error("command line argument `--" + option + "' not recognised");
          }
        }
        else
        {
          std::string const& long_option = (d.m_options.find(option))->first;

          interface_description::option_descriptor const& descriptor =
            (d.m_options.find(long_option))->second;

          if (argument.size() == option.size() + 2)   // no argument
          {
            if (descriptor.needs_argument())
            {
              error("expected argument to option `--" + option + "'!");
            }
            else if (descriptor.m_argument.get() == nullptr)
            {
              m_options.insert(std::make_pair(long_option, ""));
            }
            else
            {
              m_options.insert(std::make_pair(long_option,
                                              descriptor.m_argument->get_default()));
            }
          }
          else
          {
            std::string option_argument(argument, option.size() + 3);

            if (!descriptor.accepts_argument())
            {
              error("did not expect argument to option `--" + option + "'");
            }
            else if (!descriptor.m_argument->validate(option_argument))
            {
              error("argument `" + option_argument + "' to option `--" + option + "' is invalid");
            }

            m_options.insert(std::make_pair(long_option, option_argument));
          }
        }

        continue;
      }
      else
      {
        try
        {
          for (std::string::size_type j = 1; j < argument.size(); ++j)
          {
            std::string option(1, argument[j]);

            // Assume that the argument is a short option
            if (d.m_short_to_long.find(argument[j]) == d.m_short_to_long.end())
            {
              error("command line argument `-" + option + "' not recognised");
            }
            else
            {
              std::string const& long_option =
                d.m_options.find(d.m_short_to_long[argument[j]])->first;

              interface_description::option_descriptor const& descriptor =
                (d.m_options.find(long_option))->second;

              if (argument.size() - j == 1)   // the last option without argument
              {
                if (!descriptor.accepts_argument())   // no argument needed
                {
                  m_options.insert(std::make_pair(long_option, ""));
                }
                else if (descriptor.needs_argument())
                {
                  if (++i != arguments.end())   // next command line argument is option argument
                  {
                    m_options.insert(std::make_pair(long_option, std::string(*i)));
                  }
                  else
                  {
                    error("expected argument to option `-" + option + "'");
                  }
                }
                else
                {
                  m_options.insert(std::make_pair(long_option, descriptor.m_argument->get_default()));
                }
              }
              else   // intermediate option or option with argument
              {
                if (d.m_options.find(long_option)->second.accepts_argument())
                {
                  std::string option_argument(argument, j + 1);

                  if (!descriptor.accepts_argument())
                  {
                    error("did not expect argument to option `-" + option + "'");
                  }
                  else if (!descriptor.m_argument->validate(option_argument))
                  {
                    error("argument `" + option_argument + "' to option `-" + option + "' is invalid");
                  }

                  // must be the last option, so take the remainder as option argument
                  m_options.insert(std::make_pair(long_option, option_argument));

                  break;
                }
                else
                {
                  m_options.insert(std::make_pair(long_option, ""));
                }
              }
            }
          }
        }
        catch (std::runtime_error&)   // parse error
        {
          if (argument.substr(1, 3) == "at-")   // does not match option pattern for aterm library
          {
            continue;
          }

          throw;
        }

        continue;
      }
    }

    m_arguments.push_back(argument);
  }
}

/// \cond INTERNAL
/**
 * Converts an array of C-style strings and a count to an STL vector of STL strings.
 *
 * \param[in] count the number of arguments
 * \param[in] arguments C-style array with arguments as C-style zero-terminated string
 **/
std::vector< std::string > command_line_parser::convert(const int count, char const* const* const arguments)
{
  std::vector< std::string > result;

  if (0 < count)
  {
    result.resize(count);

    std::vector< std::string >::reverse_iterator j = result.rbegin();

    for (char const* const* i = &arguments[count - 1]; i != &arguments[0]; --i)
    {
      *(j++) = std::string(*i);
    }
  }

  return result;
}

#ifndef __CYGWIN__ // wstring is not available for cygwin
/**
 * Converts an array of C-style strings and a count to an STL vector of STL strings.
 *
 * \param[in] count the number of arguments
 * \param[in] arguments C-style array with arguments as C-style zero-terminated string
 * \pre arguments uses UTF-8 encoding
 **/
std::vector< std::string > command_line_parser::convert(const int count, wchar_t const* const* const arguments)
{
  std::vector< std::string > result;

  if (0 < count)
  {
    std::ostringstream converter;

    result.resize(count);

    std::vector< std::string >::reverse_iterator j = result.rbegin();

    for (wchar_t const* const* i = &arguments[count - 1]; i != &arguments[0]; --i)
    {
      std::wstring argument(*i);

      *(j++) = std::string(argument.begin(), argument.end());
    }
  }

  return result;
}
/// \endcond
#endif // __CYGWIN__

/**
 * Parses a string as if it were an unparsed command line and stores it as
 * an array of C-style strings and a count to an STL vector of STL strings.
 *
 * It parses the following grammar:
 *
 * \verbatim
 *  command     ::= name white-space* argument
 *  argument    ::= ( non-white-space | block ) +
 *  block       ::= "'" ( character ) * "'" | '"' ( character )* '"'
 *  white-space ::= ' ' | '\t'
 * \endverbatim
 *
 * Where `non-white-space' is a character that is neither space nor tab;
 * and `name' is a path to a file in the local filesystem and operating
 * system path syntax.
 *
 * The strings that represent 'name' and 'argument' are identified from
 * left-to-right and stored in the output vector.
 *
 * \note This method is includes specifically for on MS Windows systems. On
 * such systems the command used to start the program is not always
 * available in another form than the unparsed command line.
 *
 * \param[in] arguments C-style zero-terminated string
 **/
std::vector< std::string > command_line_parser::parse_command_line(char const* const arguments)
{
  std::vector< std::string > result;

  if (arguments != nullptr)
  {
    char const* current = arguments;

    while (*current != '\0')
    {
      // skip initial white space
      while (*current == '\0' || *current == ' ')
      {
        ++current;
      }

      if (*current != '\0')
      {
        char const* current_argument = current;

        do
        {
          if (*current == '\'')
          {
            do
            {
              ++current;
            }
            while (*current != '\0' && *current != '\'');

            if (*current == '\'')
            {
              ++current;
            }
          }
          else if (*current == '\"')
          {
            do
            {
              ++current;
            }
            while (*current != '\0' && *current != '\"');

            if (*current == '\"')
            {
              ++current;
            }
          }
          else
          {
            ++current;
          }
        }
        while (*current != '\0' && *current != ' ');

        result.emplace_back(current_argument, current - current_argument);
      }
    }
  }

  return result;
}

std::string const& command_line_parser::option_argument(std::string const& long_identifier) const
{
  if (options.count(long_identifier) == 0)
  {
    interface_description::option_descriptor const& option(m_interface.find_option(long_identifier));

    if (option.needs_argument() || option.argument().has_default())   // mandatory argument with default value
    {
      return option.get_default();
    }
    else
    {
      throw std::logic_error("Fatal error: argument requested of unspecified option!");
    }
  }
  else if (!m_interface.m_options.find(long_identifier)->second.accepts_argument())
  {
    throw std::logic_error("Fatal error: argument requested of option that does not take an argument!");
  }

  return options.find(long_identifier)->second;
}

/**
 * \param[in] d the interface description
 **/
void command_line_parser::process_default_options(interface_description& d)
{
  if (d.m_options.find("cli-testing-no-duplicate-option-checking") == d.m_options.end())
  {
    for (option_map::const_iterator i = m_options.begin(); i != m_options.end(); ++i)
    {
      if (1 < m_options.count(i->first))
      {
        error("option -" + (d.long_to_short(i->first) != '\0' ?
              std::string(1, d.long_to_short(i->first)).append(", --") : "-") + i->first + " specified more than once");
      }
    }
  }

  m_continue = false;

  if (m_options.count("help"))
  {
    std::cout << d.textual_description(false);
  }
  else if (m_options.count("help-all"))
  {
    std::cout << d.textual_description(true);
  }
  else if (m_options.count("version"))
  {
    std::cout << d.version_information();
  }
  else if (m_options.count("generate-man-page"))
  {
    std::cout << d.man_page();
  }
  else if (m_options.count("generate-xml"))
  {
    d.xml_page(std::cout);
  }
  else
  {
    using action_list = std::vector<bool (*)(command_line_parser&)>;

    action_list& actions(get_registered_actions());

    m_continue = true;

    for (action_list::const_iterator i = actions.begin(); m_continue && i != actions.end(); ++i)
    {
      m_continue &= (*i)(*this);
    }
  }
}
} // namespace mcrl2::utilities

