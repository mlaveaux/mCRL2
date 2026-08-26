// Author(s): Maurice Laveaux
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pbes/evidence_lts.h
/// \brief Writes the evidence (witness/counterexample) of a solved structure graph directly as an
/// LTS, instead of the LPS reconstructed by solve_structure_graph_with_counter_example. This
/// avoids the round trip of writing an evidence LPS and then re-exploring its full state space
/// with lps2lts, since the structure graph already contains the evidence transitions explicitly.
///
/// Kept separate from solve_structure_graph.h because it needs mcrl2/lts/lts_builder.h, which
/// pulls in mcrl2/lps/explorer.h.

#ifndef MCRL2_PBES_EVIDENCE_LTS_H
#define MCRL2_PBES_EVIDENCE_LTS_H

#include "mcrl2/lts/lts_builder.h"
#include "mcrl2/pbes/solve_structure_graph.h"

namespace mcrl2::pbes_system {

class lps_to_lts_solve_structure_graph_algorithm: public solve_structure_graph_algorithm
{
  protected:
    // Determines the initial state (values for the LPS process parameters), the same way
    // lps::explorer does.
    static lps::state initial_state(const lps::specification& lpsspec, const pbes& p, data::rewrite_strategy strategy)
    {
      data::rewriter rewr(p.data(), strategy);
      data::mutable_indexed_substitution<> sigma;
      data::data_expression_vector values;
      for (const data::data_expression& e: lpsspec.initial_process().expressions())
      {
        values.push_back(rewr(e, sigma));
      }
      return lps::state(values.begin(), values.size());
    }

    // Writes the evidence encoded by V as an LTS with the given format directly to filename,
    // using the general lts_builder abstraction (the same one lps2lts uses) instead of a bespoke
    // per-format writer.
    static void create_evidence_lts(structure_graph& G,
                                     const std::set<structure_graph::index_type>& V,
                                     const lps::specification& lpsspec,
                                     const pbes& p,
                                     const pbes_equation_index& p_index,
                                     data::rewrite_strategy strategy,
                                     lts::lts_type format,
                                     const std::string& filename)
    {
      try
      {
        lps::explorer_options options;
        options.save_at_end = true;
        std::unique_ptr<lts::lts_builder> builder = lts::create_lts_builder(lpsspec, options, format, filename);

        // Every transition encoded by V is reachable from the initial state by construction: G is
        // the winning region of the initial vertex, so it cannot contain anything else. Hence the
        // transitions can be fed to the builder directly, one by one, as long as the initial state
        // is the first one inserted: lts_builder::finalize hardcodes set_initial_state(0), so index
        // 0 must be the initial state.
        lts::lts_builder::indexed_set_for_states_type state_map;
        auto [init_index, init_inserted] = state_map.insert(initial_state(lpsspec, p, strategy));
        (void)init_index;
        (void)init_inserted;

        detail::for_each_evidence_transition(G, V, lpsspec, p, p_index,
          [&](std::size_t /* summand_index */,
              const data::data_expression_vector& source,
              const lps::multi_action& a,
              const data::data_expression_vector& target)
          {
            auto [from, from_inserted] = state_map.insert(lps::state(source.begin(), source.size()));
            auto [to, to_inserted] = state_map.insert(lps::state(target.begin(), target.size()));
            (void)from_inserted;
            (void)to_inserted;
            builder->add_transition(from, a, to);
          });

        builder->finalize(state_map, false);
        builder->save(filename);
      }
      catch (const std::exception& e)
      {
        throw mcrl2::runtime_error(std::string("Counter-example cannot be reconstructed, either wrong file provided or an internal error occurred. ") + e.what());
      }
    }

  public:
    lps_to_lts_solve_structure_graph_algorithm() = default;

    /// \brief Solve a pbes for some equation, directly writing a counter example or witness for
    /// it as an LTS, based on the accompanying linear process.
    /// \param G       A structure graph.
    /// \param lpsspec The original LPS that was used to create the PBES.
    /// \param p       The pbes to be solved.
    /// \param p_index The index of the pbes equation to be solved.
    /// \param strategy The rewrite strategy used to compute the initial state and (indirectly, via
    ///        pbesinst) the states occurring in the structure graph; must match the strategy used
    ///        to build G, so that both sides normalize data expressions to the same term.
    /// \param format  The LTS format to write \a filename in.
    /// \param filename The file to write the evidence LTS to.
    /// \return A boolean indicating the solution.
    bool solve_with_evidence_lts(structure_graph& G,
                                  const lps::specification& lpsspec,
                                  const pbes& p,
                                  const pbes_equation_index& p_index,
                                  data::rewrite_strategy strategy,
                                  lts::lts_type format,
                                  const std::string& filename)
    {
      if (!lpsspec.global_variables().empty())
      {
        throw mcrl2::runtime_error("solve_with_evidence_lts requires an LPS without global variables.");
      }
      if (!p.global_variables().empty())
      {
        throw mcrl2::runtime_error("solve_with_evidence_lts requires a PBES without global variables.");
      }
      if (lpsspec.process().has_time())
      {
        throw mcrl2::runtime_error("Generating an evidence LTS directly is not supported for timed linear processes.");
      }

      auto [is_disjunctive, W] = extract_evidence(G);
      create_evidence_lts(G, W, lpsspec, p, p_index, strategy, format, filename);
      return is_disjunctive;
    }
};

/// \brief Solve this pbes_system using a structure graph, directly writing the evidence as an LTS.
/// \param G       The structure graph.
/// \param lpsspec The original LPS that was used to create the PBES.
/// \param p       The pbes to be solved.
/// \param p_index The index of the pbes equation to be solved.
/// \param strategy The rewrite strategy used to solve \a G.
/// \param format  The LTS format to write \a filename in.
/// \param filename The file to write the evidence LTS to.
inline
bool solve_structure_graph_with_evidence_lts(structure_graph& G,
                                              const lps::specification& lpsspec,
                                              const pbes& p,
                                              const pbes_equation_index& p_index,
                                              data::rewrite_strategy strategy,
                                              lts::lts_type format,
                                              const std::string& filename)
{
  lps_to_lts_solve_structure_graph_algorithm algorithm;
  return algorithm.solve_with_evidence_lts(G, lpsspec, p, p_index, strategy, format, filename);
}

} // namespace mcrl2::pbes_system

#endif // MCRL2_PBES_EVIDENCE_LTS_H
