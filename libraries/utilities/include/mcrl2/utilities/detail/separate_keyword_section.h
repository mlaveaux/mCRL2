// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/utilities/detail/separate_keyword_section.h
/// \brief add your file description here.

#ifndef MCRL2_UTILITIES_DETAIL_SEPARATE_KEYWORD_SECTION_H
#define MCRL2_UTILITIES_DETAIL_SEPARATE_KEYWORD_SECTION_H

#include <utility>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/trim.hpp>

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/text_utility.cxx"
#else
  import utilities;
#endif

namespace mcrl2::utilities::detail
{

// Separates all sections with a keyword from the other keyword sections
// Returns a pair containing consisiting of the keyword section and the other keyword sections
inline std::pair<std::string, std::string> separate_keyword_section(const std::string& text1,
    const std::string& keyword,
    const std::vector<std::string>& all_keywords,
    bool repeat_keyword = false)
{
  std::string text = boost::trim_copy(text1);
  std::ostringstream out1; // will contain the keyword sections
  std::ostringstream out2; // will contain the other keyword sections

  std::string regex_keyword = "\\b" + keyword + "\\b";
  // create a regex that looks like this: "(\\beqn\\b)|(\\bcons\\b)|(\\bmap\\b)|(\\bvar\\b)"
  std::vector<std::string> v = all_keywords;
  v.erase(std::remove(v.begin(), v.end(), keyword), v.end()); // erase keyword from v
  for (std::string& s : v)
  {
    s = "(\\b" + s + "\\b)";
  }
  std::string regex_other_keywords = boost::algorithm::join(v, "|");

  std::vector<std::string> specs = utilities::regex_split(text, regex_keyword);
  if (!text.starts_with(keyword) && !specs.empty())
  {
    out2 << specs.front() << std::endl;
    specs.erase(specs.begin());
  }
  for (const std::string& spec : specs)
  {
    // strip trailing map/cons/var/eqn declarations
    std::vector<std::string> v = utilities::regex_split(spec, regex_other_keywords);
    if (!v.empty())
    {
      if (repeat_keyword)
      {
        out1 << "\n" << keyword;
      }
      out1 << "  " << v.front();
      out2 << spec.substr(v.front().size());
    }
  }
  std::string s1 = out1.str();
  if (!s1.empty())
  {
    if (!repeat_keyword)
    {
      s1 = keyword + "\n" + s1;
    }
  }
  return std::make_pair(s1 + "\n", out2.str() + "\n");
  }

  } // namespace mcrl2::utilities::detail

#endif // MCRL2_UTILITIES_DETAIL_SEPARATE_KEYWORD_SECTION_H
