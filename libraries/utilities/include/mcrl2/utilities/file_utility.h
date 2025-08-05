// Author(s): Sjoerd Cranen
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/utilities/file_utility.h
/// \brief add your file description here.

#ifndef MCRL2_UTILITIES_FILE_UTILITY_H
#define MCRL2_UTILITIES_FILE_UTILITY_H

#include <cassert>
#include <string>
#include <vector>

namespace mcrl2::utilities
{

/// \brief Returns true if the given filename has the extension ext. The extension does not start with a dot. 
inline
bool has_extension(const std::string& filename, const std::string& extension)
{
  assert(extension.size()>=1 && extension[0]!='.');  // The extension should not start with a dot and consist out of at least one character.
  std::string dotted_extension="."+extension;
  if (filename.size()<dotted_extension.size())
  {
    return false;
  }
  const std::string filename_extension = filename.substr(filename.size()-dotted_extension.size());
  return filename_extension==dotted_extension;
} 

/// Returns the file extension for the given path.
inline std::string file_extension(const std::string& filename)
{
  auto pos = filename.find_last_of('.');
  if (pos == std::string::npos)
  {
    return "";
  }
  return filename.substr(pos + 1);
}

class file_format
{
  protected:
    std::string m_shortname;
    std::string m_description;
    bool m_text_format;
    std::vector<std::string> m_extensions;
  
  public:
    file_format()
     : file_format("unknown", "Unknown format", false)
    {
    }

    file_format(const std::string& shortname, const std::string& description, bool is_text_format)
      : m_shortname(shortname), m_description(description), m_text_format(is_text_format)
    { }
  
    void add_extension(const std::string& ext)
    {
      m_extensions.push_back(ext);
    }
  
    const std::string& shortname() const { return m_shortname; }
  
    const std::string& description() const { return m_description; }
   
    bool text_format() const { return m_text_format; }
  
    bool matches(const std::string& filename) const
    {
      for (const std::string& ext: m_extensions)
      {
        if (has_extension(filename,ext))
        {
          return true;
        }
      }
      return false;
    }
  
    bool operator==(const file_format& other) const
    {
      return m_shortname == other.m_shortname;
    }

    bool operator<(const file_format& other) const
    {
      return m_shortname < other.m_shortname;
    }
};

inline
std::ostream& operator<<(std::ostream& stream, const file_format& format)
{
  return stream << format.shortname();
}

inline
bool file_exists(const std::string& filename)
{
  if (FILE * file = fopen(filename.c_str(), "r"))
  {
    fclose(file);
    return true;
  }
  return false;
}

} // namespace mcrl2::utilities

#endif // MCRL2_UTILITIES_FILE_UTILITY_H
