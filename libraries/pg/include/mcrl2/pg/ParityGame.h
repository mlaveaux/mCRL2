// Copyright (c) 2009-2013 University of Twente
// Copyright (c) 2009-2013 Michael Weber <michaelw@cs.utwente.nl>
// Copyright (c) 2009-2013 Maks Verver <maksverver@geocities.com>
// Copyright (c) 2009-2013 Eindhoven University of Technology
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef MCRL2_PG_PARITY_GAME_H
#define MCRL2_PG_PARITY_GAME_H

#include "mcrl2/pg/Graph.h"

#ifndef MCRL2_ENABLE_MODULES
  #include "mcrl2/utilities/exception.cxx"
#else
  import utilities;
#endif

// Forward declaration of mcrl2::pbes_system::pbes, which may or may not be
// defined later depending on whether mCRL2 support is compiled in.
 namespace mcrl2::pbes_system { class pbes; } 

#if __GNUC__ >= 3
#   define ATTR_PACKED  __attribute__((__packed__))
#else
#   define ATTR_PACKED
#endif

/*! \defgroup ParityGameData

    Data structures used to represent parity games in memory.
*/

/*! Type of priorities in the game */
 using priority_t = std::size_t;

 /*! The two players in a parity game (Even and Odd) */
 enum player_t
 {
   PLAYER_EVEN = 0, //!< Even (0)
   PLAYER_ODD = 1   //!< Odd (1)
              } ATTR_PACKED;

inline
player_t opponent(const player_t p)
{
    switch(p)
    {
      case PLAYER_EVEN:
        return PLAYER_ODD;
      case PLAYER_ODD:
        return PLAYER_EVEN;
      default:
        throw mcrl2::runtime_error("unknown player");
    }
}

/*! \ingroup ParityGameData

    Parity Game data associated with each vertex in the game graph.

    \sa ParityGame::Player
*/
struct ParityGameVertex
{
    //! the vertex owner (i.e. next player to move)
    player_t   player;

    //! the priority of the vertex  between 0 and `d` (exclusive).
    priority_t priority;
};

inline bool operator== (const ParityGameVertex &a, const ParityGameVertex &b)
{
    return a.player == b.player && a.priority == b.priority;
}

inline bool operator!= (const ParityGameVertex &a, const ParityGameVertex &b)
{
    return a.player != b.player || a.priority != b.priority;
}


/*! \ingroup ParityGameData

    A parity game extends a directed graph by assigning a player
    (Even or Odd) and an integer priority to every vertex.
    Priorities are between 0 and `d` (exclusive). */
class ParityGame
{
public:
  using Player = player_t;

  /*! A strategy determines the partitioning of the game's vertices into
      winning sets for both players and provides a deterministic strategy for
      the vertices controlled by a player in its winning set.

      For each vertex v owned by player p:
      - strategy[v] == NO_VERTEX if vertex v is not in p's winning set
      - strategy[v] == w if vertex v is in p's winning set, (v,w) is an edge
        in the game graph, and (v,w) is a winning move for player p. */
  using Strategy = std::vector<verti>;

  /*! Construct an empty parity game */
  ParityGame() = default;

  /*! Destroy a parity game */
  ~ParityGame();

  /*! Reset to an empty game. */
  void clear();

  /*! Reset the game to a copy of `game`. */
  void assign(const ParityGame& game);

  /*! Reset the game with the given graph and node attributes. */
  void assign(const StaticGraph& g, ParityGameVertex* na);

  /*! Returns whether the game is empty. */
  bool empty() const { return graph().empty(); }

  /*! Efficiently swaps the contents of this parity game with another one. */
  void swap(ParityGame& pg) noexcept;

  //!\name Generation
  //!@{

  /*! Generate a random parity game, with vertices assigned uniformly at
      random to players, and priority assigned uniformly between 0 and d-1.

      The generated game is a clustered random game if clustersize > 0, or
      an unclustered random game otherwise.

      \param V            number of game vertices
      \param clustersize  cluster size (or 0 for no clustering)
      \param outdeg       average outdegree (at least 1)
      \param edge_dir     part of edges to store
      \param d            number of priorities (at least 1)
      \sa void StaticGraph::make_random()
  */
  void make_random(verti V, unsigned clustersize, unsigned outdeg, StaticGraph::EdgeDirection edge_dir, int d);

  /*! Create a subgame containing only the given vertices from the original
      game. Vertices are renumbered to be in range [0..num_vertices).
      Edges going out of the vertex subset specified by `vertices` are
      removed, so every vertex must have at least one outgoing edge that stays
      within the vertex subset, or the result is not a valid parity game.

      \sa make_subgame(const ParityGame &, const verti *, verti, const Strategy &)
  */
  template <class ForwardIterator>
  void make_subgame(const ParityGame& game,
      ForwardIterator vertices_begin,
      ForwardIterator vertices_end,
      bool proper,
      StaticGraph::EdgeDirection edge_dir = StaticGraph::EDGE_NONE);

#ifdef MCRL2_ENABLE_MULTITHREADING
    void make_subgame_threads( const ParityGame &game,
                               const verti *verts,
                               const verti nvert,
                               bool proper,
                               StaticGraph::EdgeDirection edge_dir
                                    = StaticGraph::EDGE_NONE );
#endif // MCRL2_ENABLE_MULTITHREADING

    //!@}

    //!\name Transformation
    //!@{

    /*! Replaces the current game by its dual game, which uses the same game
        graph, but swaps players and changes priorities, such that the solution
        to the game is the same except with winners reversed. (A node won by
        even in the old game, is won by odd in the dual game, and vice versa.)
    */
    void make_dual();

    /*! Shuffles the nodes in the game. Node n becomes node perm[n], for all n.
        Requires that perm contains a permuation of 0 through graph().V()-1.

        NOTE: this is the reverse of what is normally understood to be a
              permutation! i.e. (2 3 1) maps 1->2, 2->3, 3->1 (instead of the
              usual interpretation of 1->3, 2->1, 3->2) */
    void shuffle(const std::vector<verti> &perm);

    /*! Compresses the range of priorities such that after compression,
        cardinality(p) &gt; 0, for 0 &lt; p &lt; d.

        If `cardinality` is 0, then the priorities for the game itself are used.
        Otherwise, `cardinality` must be an array of length `d` and the caller
        must ensure that all priorities that occur in the game have a positive
        cardinality count.

        If `preserve_parity` is true, then remapping preserves the parity of
        priorities and thus players of vertices. In this case, cardinality(0)
        may be zero afterwards.

        If `preserve_parity` is false, then players as well as priorities are
        remapped to preserve winning sets. The function returns the parity of
        the priority that was mapped to zero.
    */
    void compress_priorities(const verti cardinality[] = nullptr, bool preserve_parity = true);

    /*! Propagate priorities in the graph, by changing the priority for a vertex
        to the maximum of the priorities of its successors, if this is less than
        the vertex's current priority value, and similarly for its predecessors.

        This preserves winners and optimal strategies, but usually shifts the
        distribution of priorities towards lower values, which benefits the
        performance in some solvers.

        Returns the sum of differences between old and new priority values
        for all vertices.
    */
    long long propagate_priorities();
    //!@}

    //!\name Input/Output
    //!@{

    /*! Read a game description in PGSolver format. */
    void read_pgsolver( std::istream &is,
        StaticGraph::EdgeDirection edge_dir = StaticGraph::EDGE_BIDIRECTIONAL );

    /*! Write a game description in PGSolver format. */
    void write_pgsolver(std::ostream &os) const;

    /*! Read a game description from an mCRL2 PBES. */
    void read_pbes(const std::string& file_path,
        verti* goal_vertex = nullptr,
        StaticGraph::EdgeDirection edge_dir = StaticGraph::EDGE_BIDIRECTIONAL,
        const std::string& rewrite_strategy = "jitty");

    /*! Read raw parity game data from input stream */
    void read_raw(std::istream &is);

    /*! Write raw parity game data to output stream */
    void write_raw(std::ostream &os) const;

    /*! Write a game description in Graphviz DOT format */
    void write_dot(std::ostream &os) const;

    /*! Write human-readable description of game to error stream (intended to
        be used while debugging only) with optional strategy. */
    void write_debug( const Strategy &s = Strategy(),
        std::ostream &os = std::cerr) const;

    //!@}

    //! \name Data access
    //!@{

    /*! Returns the priority limit d; all priorities must be in range [0:d). */
    priority_t d() const { return d_; }

    /*! Returns the game graph */
    const StaticGraph &graph() const { return graph_; }

    /*! Returns the priority associated with vertex v */
    priority_t priority(verti v) const { return vertex_[v].priority; }

    /*! Returns the player associated with vertex v */
    Player player(verti v) const { return (Player)(signed char)vertex_[v].player; }

    /*! Returns the number of vertices with priority `p`.
        `p` must be between 0 and `d` (exclusive). */
    verti cardinality(int p) const { return cardinality_[p]; }

    /*! Returns the winner for vertex v according to strategy s. */
    template<class StrategyT>
    Player winner(const StrategyT &s, verti v) const;

    //!@}

    //!\name Verification
    //!@{

    /*! Checks if this is a proper game.

        A game is proper if every vertex has a successor in the game graph. */
    bool proper() const;

    /*! Verifies that the given strategy is valid in this game.

        \param s the strategy to be verified
        \param error if non-NULL, then *error is set NO_VERTEX if the strategy is
                     correct, and to the index of an incorrectly classified
                     vertex if the strategy is incorrect.
    */
    bool verify(const Strategy &s, verti *error) const;
    //!@}

    /*! Generate a parity game from an mCRL2 PBES. */
    void assign_pbes(mcrl2::pbes_system::pbes& pbes,
        verti* goal_vertex = nullptr,
        StaticGraph::EdgeDirection edge_dir = StaticGraph::EDGE_BIDIRECTIONAL,
        const std::string& rewrite_strategy = "jitty");

  protected:
    /*! Re-allocate memory to store information on V vertices with priorities
        between 0 and `d` (exclusive). */
    void reset(verti V, int d);

    /*! Recalculate cardinalities (priority frequencies) from the first
        `num_vertices` elements of `vertex_`. */
    void recalculate_cardinalities(verti num_vertices);

    /*! Helper function for ParityGame::propagate_priorities() that decreases
        the priority for `v` to the maximum of those in range [begin:end), if
        this is less than its current value, and returns the absolute change. */
    int propagate_priority( verti v, StaticGraph::const_iterator it,
                                     StaticGraph::const_iterator end );

private:
    explicit ParityGame(const ParityGame &game);
    ParityGame &operator=(const ParityGame &game);

private:
  int d_ = 0;         /*!< priority limit (max. priority + 1) */
  StaticGraph graph_; /*!< game graph */

  /*! Assignment of players and priorities to vertices (size graph_.V()) */
  ParityGameVertex* vertex_ = nullptr;

  /*! Cardinality counts for priorities.
      cardinality_[p] is equal to the number of vertices with priority p. */
  verti* cardinality_ = nullptr;
};

inline void swap(ParityGame& a, ParityGame& b) noexcept
{
  a.swap(b);
}

#include "ParityGame_impl.h"

#endif /* ndef MCRL2_PG_PARITY_GAME_H */
