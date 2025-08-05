// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pbes/pbespgsolve.h
/// \brief add your file description here.

#ifndef MCRL2_PBES_PBESPGSOLVE_H
#define MCRL2_PBES_PBESPGSOLVE_H

#include <memory>

#include "mcrl2/data/rewriter.h"
#include "mcrl2/pbes/algorithms.h"
#include "mcrl2/pg/ComponentSolver.h"
#include "mcrl2/pg/DecycleSolver.h"
#include "mcrl2/pg/DeloopSolver.h"
#include "mcrl2/pg/PredecessorLiftingStrategy.h"
#include "mcrl2/pg/PriorityPromotionSolver.h"

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/execution_timer.cxx"
#else
  import utilities;
#endif

namespace mcrl2::pbes_system {

enum pbespg_solver_type
{
  spm_solver,
  alternative_spm_solver,
  recursive_solver,
  priority_promotion
};

inline
pbespg_solver_type parse_solver_type(const std::string& s)
{
  if (s == "spm")
  {
    return spm_solver;
  }
  else if (s == "altspm")
  {
    return alternative_spm_solver;
  }
  else if (s == "recursive")
  {
    return recursive_solver;
  }
  else if (s == "prioprom")
  {
    return priority_promotion;
  }
  throw mcrl2::runtime_error("unknown solver " + s);
}

/// \return A string representation of the solver type.
inline
std::string print(const pbespg_solver_type solver_type)
{
  switch(solver_type)
  {
    case spm_solver: return "spm";
    case alternative_spm_solver: return "altspm";
    case recursive_solver: return "recursive";
    case priority_promotion: return "prioprom";
  }
  throw mcrl2::runtime_error("unknown solver");
}

// \return A string representation of the solver type.
inline
std::string description(const pbespg_solver_type solver_type)
{
  switch(solver_type)
  {
    case spm_solver: return "Small progress measures";
    case alternative_spm_solver: return "Alternative implementation of small progress measures";
    case recursive_solver: return "Recursive algorithm";
    case priority_promotion: return "Priority promotion (experimental)";
  }
  throw mcrl2::runtime_error("unknown solver");
}

inline
std::istream& operator>>(std::istream& is, pbespg_solver_type& t)
{
  try
  {
    std::string s;
    is >> s;
    t = parse_solver_type(s);
  }
  catch(mcrl2::runtime_error&)
  {
    is.setstate(std::ios_base::failbit);
  }
  return is;
}

inline
std::ostream& operator<<(std::ostream& os, const pbespg_solver_type t)
{
  os << print(t);
  return os;
}


/// \return A string representation of the edge direction.
inline
std::string print(StaticGraph::EdgeDirection edge_direction)
{
  if (edge_direction == StaticGraph::EDGE_SUCCESSOR)
  {
    return "successor";
  }
  else if (edge_direction == StaticGraph::EDGE_PREDECESSOR)
  {
    return "predecessor";
  }
  else if (edge_direction == StaticGraph::EDGE_BIDIRECTIONAL)
  {
    return "bidirectional";
  }
  return "unknown edge direction";
}

struct pbespgsolve_options
{
  pbespg_solver_type solver_type = spm_solver;
  bool use_scc_decomposition = true;
  bool use_decycle_solver = false;
  bool use_deloop_solver = true;
  bool verify_solution = true;
  bool only_generate = false;
  data::rewriter::strategy rewrite_strategy = data::jitty;
};

class pbespgsolve_algorithm
{
  protected:
    std::unique_ptr<ParityGameSolverFactory> solver_factory;
    mcrl2::utilities::execution_timer& m_timer;
    pbespgsolve_options m_options;

  public:

    pbespgsolve_algorithm(mcrl2::utilities::execution_timer& timing, pbespgsolve_options options = pbespgsolve_options())
      : m_timer(timing),
        m_options(options)
    {
      if (options.solver_type == spm_solver || options.solver_type == alternative_spm_solver)
      {
        bool alternative_solver = (options.solver_type == alternative_spm_solver);

        // Create a SPM solver factory:
        solver_factory = std::make_unique<SmallProgressMeasuresSolverFactory>(
            std::make_shared<PredecessorLiftingStrategyFactory>(),
            2,
            alternative_solver);
      }
      else if (options.solver_type == recursive_solver)
      {
        // Create a recursive solver factory:
        solver_factory = std::make_unique<RecursiveSolverFactory>();
      }
      else if (options.solver_type == priority_promotion)
      {
        solver_factory = std::make_unique<PriorityPromotionSolverFactory>();
      }
      else
      {
        throw mcrl2::runtime_error("pbespgsolve: unknown solver type");
      }

      if (options.use_scc_decomposition)
      {
        // Wrap solver factory into a component solver factory:
        solver_factory = std::make_unique<ComponentSolverFactory>(*solver_factory.release());
      }

      if (options.use_decycle_solver)
      {
        solver_factory = std::make_unique<DecycleSolverFactory>(*solver_factory.release());
      }

      if (options.use_deloop_solver)
      {
        solver_factory = std::make_unique<DeloopSolverFactory>(*solver_factory.release());
      }
    }

    bool run(ParityGame& pg, const verti goal_v)
    {
      if (!m_options.only_generate)
      {
        mCRL2log(log::verbose) << "Solving..." << std::endl;

        // Create a solver:
        std::unique_ptr<ParityGameSolver> solver(solver_factory->create(pg));

        // Solve the game:
        m_timer.start("solving");
        ParityGame::Strategy solution = solver->solve();
        m_timer.finish("solving");

        if (solution.empty())
        {
          throw mcrl2::runtime_error("pbespgsolve: solving failed!\n");
        }

        verti error_vertex;

        // Optional: verify the solution
        if (m_options.verify_solution && !pg.verify(solution, &error_vertex))
        {
          throw mcrl2::runtime_error("pbespgsolve: verification of the solution failed!\n");
        }

        return pg.winner(solution, goal_v) == PLAYER_EVEN;
      }
      return true;
    }

    bool run(pbes& p)
    {
      m_timer.start("initialization");
      mCRL2log(log::verbose) << "Generating parity game..."  << std::endl;
      // Generate the game from a PBES:
      verti goal_v;
      ParityGame pg;

      pg.assign_pbes(p, &goal_v, StaticGraph::EDGE_BIDIRECTIONAL, data::pp(m_options.rewrite_strategy)); // N.B. mCRL2 could raise an exception here
      mCRL2log(log::verbose) << "Game: " << pg.graph().V() << " vertices, " << pg.graph().E() << " edges." << std::endl;
      m_timer.finish("initialization");

      return run(pg, goal_v);
    }
};


/// \brief Solves a pbes using a parity game solver
/// \return The solution of the pbes
inline
bool pbespgsolve(pbes& p, const pbespgsolve_options& options = pbespgsolve_options())
{
  utilities::execution_timer timer;
  pbes_system::algorithms::instantiate_global_variables(p);
  pbespgsolve_algorithm algorithm(timer, options);
  return algorithm.run(p);
}

/// \brief Solves a pbes using a parity game solver
/// \return The solution of the pbes
inline
bool pbespgsolve(pbes& p, utilities::execution_timer& timer, const pbespgsolve_options& options = pbespgsolve_options())
{
  pbes_system::algorithms::instantiate_global_variables(p);
  pbespgsolve_algorithm algorithm(timer, options);
  return algorithm.run(p);
}

} // namespace mcrl2::pbes_system



#endif // MCRL2_PBES_PBESPGSOLVE_H
