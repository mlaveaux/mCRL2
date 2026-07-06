// Author(s): Luc Engelen
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/data/detail/prover/bdd2dot.h
/// \brief Interface to class BDD2Dot

#ifndef MCRL2_DATA_DETAIL_PROVER_BDD2DOT_H
#define MCRL2_DATA_DETAIL_PROVER_BDD2DOT_H

#include "mcrl2/data/detail/prover/bdd_info.h"
#include <fstream>

namespace mcrl2::data::detail
{

/// \brief The class BDD2Dot offers the ability to write binary decision diagrams to dot files.
/// The method BDD2Dot::output_bdd receives a binary decision diagram as parameter a_bdd and writes it to a file in dot
/// format with the name passed as parameter a_file_name.

class BDD2Dot
{
  private:
    /// \brief The smallest unused node number.
    int f_node_number = 0;

    /// \brief The file the output is written to.
    std::ofstream f_dot_file;

    /// \brief A table containing all the visited nodes. It maps these nodes to the corresponding node numbers.
    std::map < atermpp::aterm, atermpp::aterm_int> f_visited;

    /// \brief A class that gives information about the structure of BDDs.
    BDD_Info f_bdd_info;

    /// \brief Writes the BDD it receives to BDD2Dot::f_dot_file.
    /// \param a_bdd A binary decision diagram.
    /// \pre The argument passed as parameter a_bdd is a data expression in internal mCRL2 format with the
    /// following restrictions: It either represents the constant true or the constant false, or it is an if-then-else
    /// expression with an expression of Bool as guard, and a then-branch and an else-branch that again follow these
    /// restrictions
    void aux_output_bdd(const data_expression &a_bdd)
    {
      if (f_visited.count(a_bdd)>0)  // a_bdd has already been visited.
      {
        return;
      }

      if (f_bdd_info.is_true(a_bdd))
      {
        f_dot_file << "  " << f_node_number << " [shape=box, label=\"T\"];" << std::endl;
      }
      else if (f_bdd_info.is_false(a_bdd))
      {
        f_dot_file << "  " << f_node_number << " [shape=box, label=\"F\"];" << std::endl;
      }
      else if (f_bdd_info.is_if_then_else(a_bdd))
      {
        const data_expression v_true_branch = f_bdd_info.get_true_branch(a_bdd);
        const data_expression v_false_branch = f_bdd_info.get_false_branch(a_bdd);
        aux_output_bdd(v_true_branch);
        aux_output_bdd(v_false_branch);
        std::size_t v_true_number = f_visited[v_true_branch].value();
        std::size_t v_false_number = f_visited[v_false_branch].value();
        const data_expression v_guard = f_bdd_info.get_guard(a_bdd);
        f_dot_file << "  " << f_node_number << " [label=\"" << mcrl2::data::pp(v_guard) << "\"];" << std::endl;
        f_dot_file << "  " << f_node_number << " -> " << v_true_number << ";" << std::endl;
        f_dot_file << "  " << f_node_number << " -> " << v_false_number << " [style=dashed];" << std::endl;
      }
      else
      {
        f_dot_file << "  " << f_node_number << " [shape=box, label=\"" << mcrl2::data::pp(a_bdd) << "\"];" << std::endl;
      }
      f_visited[a_bdd]= atermpp::aterm_int(f_node_number++);
    }

  public:
    /// \brief Writes the BDD it receives as input to a file
    /// \brief in dot format with the name received as input.
    /// precondition: The argument passed as parameter a_bdd is a data expression in internal mCRL2 format with the
    /// following restrictions: It either represents the constant true or the constant false, or it is an if-then-else
    /// expression with an expression of Bool as guard, and a then-branch and an else-branch that again follow these
    /// restrictions
    /// \param a_bdd A binary decision diagram.
    /// \param a_file_name A file name.
    void output_bdd(const data_expression &a_bdd, const std::string& a_file_name)
    {
      f_node_number = 0;
      f_dot_file.open(a_file_name);
      f_dot_file << "digraph BDD {" << std::endl;
      aux_output_bdd(a_bdd);
      f_dot_file << "}" << std::endl;
      f_dot_file.close();
    }
};

} // namespace mcrl2::data::detail

#endif
