// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include "mcrl2/data/detail/rewrite/innermost.h"

#include "mcrl2/utilities/logger.h"

#include <assert.h>

using namespace mcrl2::data;
using namespace mcrl2::data::detail;

using namespace mcrl2::log;

InnermostRewriter::InnermostRewriter(const data_specification& data_spec, const used_data_equation_selector& selector)
  : Rewriter(data_spec, selector)
{
  // This probably belongs inside the Rewriter class, checks for all the selected equations whether they are valid rewrite rules.
  for (const data_equation& eq: data_spec.equations())
  {
    if (selector(eq))
    {
      try
      {
        CheckRewriteRule(eq);
      }
      catch (std::runtime_error& e)
      {
        mCRL2log(warning) << e.what() << std::endl;
        continue;
      }
    }
  }
}

data_expression InnermostRewriter::rewrite(const data_expression& term, substitution_type& sigma)
{
  return rewrite_impl(term, sigma, {});
}

// Private functions

data_expression InnermostRewriter::rewrite_impl(const data_expression& term, substitution_type& sigma, std::set<variable> free_variables)
{
  // If t in variables
  if (is_variable(term))
  {
    const auto& var = static_cast<const data::variable&>(term);
    return sigma(var);
  }
  // Else if t in function_symbols
  else if (is_function_symbol(term))
  {
    return term;
  }
  // Else if t is of the form lambda x . u
  else if (is_abstraction(term))
  {
    const auto& abstraction = static_cast<const data::abstraction&>(term);
    assert(false);
    return term;
  }
  // Else (t is of the form f(u_0, ..., u_n)).
  else
  {

  }

  assert(false);
  return term;
}
