// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/enumerator_identifier_generator.h
/// \brief The class enumerator_identifier_genierator 
///        contains a generator of fresh identifier strings that
///        works intimately together with the atermpp function_symbols
///        to make it as efficient to generate fresh strings. Construction of a
///        enumerator identifier generator is less efficient, as it sets up
///        internal data structures and searches through all existing
///        function symbols in the aterm library. 

#ifndef MCRL2_DATA_DETAIL_ENUMERATOR_IDENTIFIER_GENERATOR_H
#define MCRL2_DATA_DETAIL_ENUMERATOR_IDENTIFIER_GENERATOR_H

#include "mcrl2/core/identifier_string.h"

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/atermpp/function_symbol_generator.cxx"
#else
  import atermpp;
#endif

namespace mcrl2::data
{

class enumerator_identifier_generator
{
  protected:
    atermpp::function_symbol_generator f;

  public:
    /// \brief Constructor
    /// \param prefix The prefix of the generated generated strings
    /// \pre The prefix may not be empty, and it may not have trailing digits
    explicit enumerator_identifier_generator(const std::string& prefix = "x_")
      : f(prefix)
    { }
    
    /// \brief Generates a unique function symbol with the given prefix followed by a number.
    core::identifier_string operator()()
    {
      return core::identifier_string(f());
    }
    
    void clear()
    {
      f.clear();
    }
};

} // namespace mcrl2::data

#endif // MCRL2_DATA_DETAIL_ENUMERATOR_IDENTIFIER_GENERATOR_H
