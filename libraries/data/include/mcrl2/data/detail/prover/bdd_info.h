// Author(s): Luc Engelen
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/detail/prover/bdd_info.h
/// \brief Interface to class BDD_Info

#ifndef MCRL2_DATA_DETAIL_PROVER_BDD_INFO_H
#define MCRL2_DATA_DETAIL_PROVER_BDD_INFO_H

#include "mcrl2/data/bool.h"

namespace mcrl2::data::detail
{

/// \brief The class BDD_Info provides information about the structure of binary decision diagrams.
class BDD_Info
{
  protected:
    static const mcrl2::data::data_expression& argument(const mcrl2::data::data_expression& x, std::size_t n)
    {
      const mcrl2::data::application& a = atermpp::down_cast<mcrl2::data::application>(x);
      return a[n]; 
    }

  public:

    /// \brief Method that returns the guard of a BDD.
    /// \param[in] a_bdd A binary decision diagram.
    /// \return The guard at the root of the BDD.
    static inline const mcrl2::data::data_expression& get_guard(const mcrl2::data::data_expression& a_bdd)
    {
      return argument(a_bdd,0);
    }

    /// \brief Method that returns the true-branch of a BDD.
    /// \param[in] a_bdd A binary decision diagram.
    /// \return The true-branch of the BDD.
    static inline const mcrl2::data::data_expression& get_true_branch(const mcrl2::data::data_expression&  a_bdd)
    {
      return argument(a_bdd,1);
    }

    /// \brief Method that returns the false-branch of a BDD.
    /// \param[in] a_bdd A binary decision diagram.
    /// \return The false-branch of the BDD.
    static inline const mcrl2::data::data_expression& get_false_branch(const mcrl2::data::data_expression& a_bdd)
    {
      return argument(a_bdd,2);
    }

    /// \brief Method that indicates whether or not a BDD equals true.
    /// \param[in] a_bdd A binary decision diagram.
    /// \return True, if the BDD equals true.
    ///         False, if the BDD does not equal true.
    static inline bool is_true(const data_expression& a_bdd)
    {
      return mcrl2::data::sort_bool::is_true_function_symbol(a_bdd);
    }

    /// \brief Method that indicates whether or not a BDD equals false.
    /// \param[in] a_bdd A binary decision diagram.
    /// \return True, if the BDD equals false.
    ///         False, if the BDD does not equal true.
    static inline bool is_false(const data_expression& a_bdd)
    {
      return mcrl2::data::sort_bool::is_false_function_symbol(a_bdd);
    }

    /// \brief Method that indicates wether or not the root of a BDD is a guard node.
    /// \param[in] a_bdd A binary decision diagram.
    /// \return True, if the root of the BDD is a guard node.
    ///         False, if the BDD equals true or if the BDD equals false.
    static inline bool is_if_then_else(const data_expression& a_bdd)
    {
      return mcrl2::data::is_if_application(a_bdd);
    }
};

} // namespace mcrl2::data::detail

#endif
