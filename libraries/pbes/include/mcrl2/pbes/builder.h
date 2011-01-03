// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://svn.win.tue.nl/trac/MCRL2/browser/trunk/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pbes/builder.h
/// \brief add your file description here.

#ifndef MCRL2_PBES_PBES_H
#include "mcrl2/pbes/pbes.h"
#endif

#ifndef MCRL2_PBES_BUILDER_H
#define MCRL2_PBES_BUILDER_H

#include "mcrl2/core/builder.h"
#include "mcrl2/data/builder.h"
#include "mcrl2/pbes/pbes.h"

namespace mcrl2 {

namespace pbes_system {

  /// \brief Traversal class for pbes expressions
  template <typename Derived>
  class pbes_expression_builder: public core::builder<Derived>
  {
    public:
      typedef core::builder<Derived> super;
      using super::enter;
      using super::leave;
      using super::operator();

#include "mcrl2/pbes/detail/pbes_expression_builder.inc.h"
  };

  /// \brief Traversal class for data expressions
  template <typename Derived>
  class data_expression_builder: public data::data_expression_builder<Derived>
  {
    public:
      typedef data::data_expression_builder<Derived> super;
      using super::enter;
      using super::leave;
      using super::operator();

#include "mcrl2/pbes/detail/data_expression_builder.inc.h"
  };

  /// \brief Traversal class for data expressions
  template <typename Derived>
  class sort_expression_builder: public data::sort_expression_builder<Derived>
  {
    public:
      typedef data::sort_expression_builder<Derived> super;
      using super::enter;
      using super::leave;
      using super::operator();

#include "mcrl2/pbes/detail/sort_expression_builder.inc.h"
  };

} // namespace pbes_system

} // namespace mcrl2

#endif // MCRL2_PBES_BUILDER_H
