// Copyright (c) 2009-2013 University of Twente
// Copyright (c) 2009-2013 Michael Weber <michaelw@cs.utwente.nl>
// Copyright (c) 2009-2013 Maks Verver <maksverver@geocities.com>
// Copyright (c) 2009-2013 Eindhoven University of Technology
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "mcrl2/pg/ParityGame.h"
#include "mcrl2/pg/SCC.h"

/*! \file pg/source/ParityGame_verify.cpp

    Implementation of parity game solution verification algorithm.
*/


/*! A functor that is used to verify the solution for a strongly-connected
    component of the game graph.

    This is a separate class because the SCC decomposition algorithm expects
    a functor which is called for every component found.

    \see ParityGame::verify
*/
struct VerifySCC  // used by ParityGame::verify()
{
    const ParityGame    &game;
    const StaticGraph   &graph;
    const priority_t    prio;
    verti * const       error;

    int operator() (const verti *scc, std::size_t scc_size)
    {
        // Search vertices in this SCC for a vertex with priority `prio`:
        for (std::size_t i = 0; i < scc_size; ++i)
        {
            verti v = scc[i];
            if (game.priority(v) == prio)
            {
                // Cycle detected if |SCC| > 1 or v has a self-edge:
                if (scc_size > 1 || graph.has_succ(v, v))
                {
                  if (error)
                  {
                    *error = v;
                  }
                    return 1;
                }
            }
        }
        return 0;
    }
};

bool ParityGame::verify(const Strategy &s, verti *error) const
{
    assert(s.size() == graph_.V());

    /* Make sure winning sets are consistently defined; i.e. only existent
       edges are used, and there are no transitions that cross winning sets. */
    for (verti v = 0; v < graph_.V(); ++v)
    {
        Player pl = winner(s, v);

        if (pl == player(v))  /* vertex won by owner */
        {
            /* Verify that:
                1. the vertex owner has a strategy (necessarily passes)
                2. the strategy uses an existing edge
                3. the strategy doesn't move outside the owner's winning set */
            if ( s[v] == NO_VERTEX || !graph_.has_succ(v, s[v]) ||
                 winner(s, s[v]) != pl )
            {
              if (error)
              {
                *error = v;
              }
                return false;
            }
        }
        else  /* vertex lost by owner */
        {
            // Verify owner has no strategy: (necessarily passes)
            if (s[v] != NO_VERTEX)
            {
              if (error)
              {
                *error = v;
              }
                return false;
            }

            // Verify owner cannot move outside this winning set:
            for (StaticGraph::const_iterator it = graph_.succ_begin(v);
                 it != graph_.succ_end(v); ++it)
            {
                if (winner(s, *it) != pl)
                {
                  if (error)
                  {
                    *error = v;
                  }
                    return false;
                }
            }
        }
    }

    // Verify absence of cycles owned by opponent in winning sets
    for (priority_t prio = 0; prio < d(); ++prio)
    {
        /* Create set of edges incident with vertices in the winning set of
           player (1 - prio%2) consistent with strategy s and incident with
           vertices of priorities >= prio only. */
        StaticGraph::edge_list edges;
        for (verti v = 0; v < graph_.V(); ++v)
        {
            if (priority(v) >= prio && (int)winner(s, v) == (1 - prio%2))
            {
                if (s[v] != NO_VERTEX)
                {
                    if (priority(s[v]) >= prio)
                    {
                      edges.emplace_back(v, s[v]);
                    }
                }
                else
                {
                    for (StaticGraph::const_iterator it = graph_.succ_begin(v);
                         it != graph_.succ_end(v); ++it)
                    {
                        if (priority(*it) >= prio)
                        {
                          edges.emplace_back(v, *it);
                        }
                    }
                }
            }
        }

        /* NOTE: we should NOT compact vertices here, because then we cannot
           use their indices to determine the priority of vertices in
           VerifySCC::operator().

           Alternatively, we could compress the priority vector as well.
           We could share some code with make_subgame() to do this (except that
           we don't need the player vector, and we don't need bidirectional
           edges, which require sorting).
        */

        // Create a subgraph storing successors only:
        StaticGraph subgraph;
        subgraph.assign(edges, StaticGraph::EDGE_SUCCESSOR);

        // Find a vertex with priority prio on a cycle:
        VerifySCC verifier = {.game = *this, .graph = subgraph, .prio = prio, .error = error};
        if (decompose_graph(subgraph, verifier) != 0)
        {
            // VerifySCC has already set *error here.
            return false;
        }
    }
    if (error)
    {
      *error = NO_VERTEX;
    }
    return true;
}
