// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pbes/pbesinst_structure_graph.h
/// \brief A variant of the lazy algorithm for instantiating a PBES, that produces a structure_graph.

#ifndef MCRL2_PBES_PBESINST_STRUCTURE_GRAPH_H
#define MCRL2_PBES_PBESINST_STRUCTURE_GRAPH_H

#include <iomanip>

#include "mcrl2/pbes/algorithms.h"
#include "mcrl2/pbes/join.h"
#include "mcrl2/pbes/pbes.h"
#include "mcrl2/pbes/pbes_quotient.h"
#include "mcrl2/pbes/pbesinst_lazy.h"
#include "mcrl2/pbes/structure_graph_builder.h"
#include "mcrl2/utilities/logger.h"



namespace mcrl2::pbes_system {

/// \brief Variant of pbesinst that will compute a structure graph for a PBES.
/// The result will be put in the structure graph that is passed in the constructor.
class pbesinst_structure_graph_algorithm: public pbesinst_lazy_algorithm
{
  protected:
    detail::structure_graph_builder m_graph_builder;

    void SG0(const propositional_variable_instantiation& X, const pbes_expression& psi, std::size_t k)
    {
      auto vertex_phi = m_graph_builder.insert_variable(X, psi, k);
      if (is_true(psi))
      {
        // skip
      }
      else if (is_false(psi))
      {
        // skip
      }
      else if (is_propositional_variable_instantiation(psi))
      {
        auto vertex_psi = m_graph_builder.insert_variable(psi);
        m_graph_builder.insert_edge(vertex_phi, vertex_psi);
      }
      else if (is_and(psi))
      {
        for (const pbes_expression& psi_i: split_and(psi))
        {
          auto vertex_psi_i = SG1(psi_i);
          m_graph_builder.insert_edge(vertex_phi, vertex_psi_i);
        }
      }
      else if (is_or(psi))
      {
        for (const pbes_expression& psi_i: split_or(psi))
        {
          auto vertex_psi_i = SG1(psi_i);
          m_graph_builder.insert_edge(vertex_phi, vertex_psi_i);
        }
      }
    }

    structure_graph::index_type SG1(const pbes_expression& psi)
    {
      auto vertex_psi = m_graph_builder.insert_vertex(psi);
      if (is_true(psi))
      {
        // skip
      }
      else if (is_false(psi))
      {
        // skip
      }
      else if (is_propositional_variable_instantiation(psi))
      {
        // skip
      }
      else if (is_and(psi))
      {
        for (const pbes_expression& psi_i: split_and(psi))
        {
          auto vertex_psi_i = SG1(psi_i);
          m_graph_builder.insert_edge(vertex_psi, vertex_psi_i);
        }
      }
      else if (is_or(psi))
      {
        for (const pbes_expression& psi_i: split_or(psi))
        {
          auto vertex_psi_i = SG1(psi_i);
          m_graph_builder.insert_edge(vertex_psi, vertex_psi_i);
        }
      }
      return vertex_psi;
    }

    std::optional<std::string> status_message(std::size_t equation_count) override
    {
      if (equation_count > 0 && equation_count % 1000 == 0)
      {
        std::ostringstream out;
        out << "Generated " << equation_count << " BES equations (" << std::fixed << std::setprecision(2) <<
          ((100.0 * equation_count) / m_graph_builder.extent()) << "% explored)" << std::endl;
        return out.str();
      }
      
      return std::nullopt;
    }

  public:
    pbesinst_structure_graph_algorithm(
      const pbessolve_options& options,
      const pbes& p,
      structure_graph& G,
      pbes_quotient& quotient,
      std::optional<data::rewriter> rewriter = std::nullopt
    )
      : pbesinst_lazy_algorithm(options, p, rewriter),
        m_graph_builder(G)
    {}

    void on_report_equation(const std::size_t /* thread_index */,
                            const propositional_variable_instantiation& X,
                            const pbes_expression& psi,
                            std::size_t k
                           ) override
    {
      // the body of this if statement will only be executed for the first equation
      if (m_graph_builder.m_initial_state == data::undefined_data_expression())
      {
        m_graph_builder.set_initial_state(X);
      }
      SG0(X, psi, k);
    }

    void run() override
    {
      pbesinst_lazy_algorithm::run();
      m_graph_builder.finalize();
      mCRL2log(log::verbose) << "Generated " << m_graph_builder.m_graph.num_of_edges() << " edges in structured graph." << std::endl;
    }
};

} // namespace mcrl2::pbes_system



#endif // MCRL2_PBES_PBESINST_STRUCTURE_GRAPH_H
