// Author(s): Sjoerd Cranen
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/core/detail/dparser_functions.h
/// \brief add your file description here.

#ifndef MCRL2_CORE_DETAIL_DPARSER_FUNCTIONS_H
#define MCRL2_CORE_DETAIL_DPARSER_FUNCTIONS_H

#include "mcrl2/core/dparser.h"

// prototypes
struct D_ParseNode;
struct D_Parser;

namespace mcrl2::core::detail
{

/// \brief Function for resolving ambiguities in the '_ -> _ <> _' operator for
/// process expressions.
struct D_ParseNode* ambiguity_fn(struct D_Parser*, int, struct D_ParseNode**);

/// \brief Custom syntax error function that prints both the line number and the column.
void syntax_error_fn(struct D_Parser *ap);

} // namespace mcrl2::core::detail

#endif // MCRL2_CORE_DETAIL_DPARSER_FUNCTIONS_H
