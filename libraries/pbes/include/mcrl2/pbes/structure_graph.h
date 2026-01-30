// Author(s): Wieger Wesselink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file mcrl2/pbes/structure_graph.h
/// \brief add your file description here.

#ifndef MCRL2_PBES_STRUCTURE_GRAPH_H
#define MCRL2_PBES_STRUCTURE_GRAPH_H

#include <iomanip>
#include <boost/dynamic_bitset.hpp>
#include <boost/range/adaptor/filtered.hpp>

#include "mcrl2/atermpp/standard_containers/vector.h"
#include "mcrl2/core/detail/print_utility.h"
#include "mcrl2/pbes/pbes.h"

namespace mcrl2::pbes_system {

namespace detail {

struct structure_graph_builder;
struct manual_structure_graph_builder;

} // namespace detail

constexpr inline
unsigned int undefined_vertex()
{
  return std::numeric_limits<unsigned int>::max();
}

// A structure graph with a facility to exclude a subset of the vertices.
// It has the same interface as simple_structure_graph.
class structure_graph
{
  friend struct detail::structure_graph_builder;
  friend struct detail::manual_structure_graph_builder;


  public:
    enum decoration_type
    {
      d_disjunction = 0,
      d_conjunction = 1,
      d_true,
      d_false,
      d_none
    };

    using index_type = unsigned int;

    struct vertex
    {
      atermpp::detail::reference_aterm<pbes_expression> m_formula;
      decoration_type decoration;
      std::size_t rank;
      std::vector<index_type> predecessors;
      std::vector<index_type> successors;
      mutable index_type strategy;

      explicit vertex(pbes_expression  formula_,
             decoration_type decoration_ = structure_graph::d_none,
             std::size_t rank_ = data::undefined_index(),
             std::vector<index_type> pred_ = std::vector<index_type>(),
             std::vector<index_type> succ_ = std::vector<index_type>(),
             index_type strategy_ = undefined_vertex()
            )
        : m_formula(std::move(formula_)),
          decoration(decoration_),
          rank(rank_),
          predecessors(std::move(pred_)),
          successors(std::move(succ_)),
          strategy(strategy_)
      {}

      void remove_predecessor(index_type u)
      {
        predecessors.erase(std::remove(predecessors.begin(), predecessors.end(), u), predecessors.end());
      }

      // Downcast reference aterm
      inline const pbes_expression& formula() const
      {
        return m_formula;
      }

      void remove_successor(index_type u)
      {
        successors.erase(std::remove(successors.begin(), successors.end(), u), successors.end());
      }

      bool is_defined() const
      {
        return  ((decoration != structure_graph::d_none) || (rank != data::undefined_index()))
             && (!successors.empty() || (decoration == d_true || decoration == d_false));
      }
      
      void inline mark(atermpp::term_mark_stack& todo) const
      {
        mark_term(m_formula, todo);
      }
    };

    using vertex_vector = atermpp::vector<vertex, std::allocator<atermpp::detail::reference_aterm<vertex>>, mcrl2::utilities::detail::GlobalThreadSafe>;

  protected:
    vertex_vector m_vertices;
    index_type m_initial_vertex = 0;
    boost::dynamic_bitset<> m_exclude;

    struct integers_not_contained_in
    {
      const boost::dynamic_bitset<>& subset;

      explicit integers_not_contained_in(const boost::dynamic_bitset<>& subset_)
        : subset(subset_)
      {}

      bool operator()(index_type i) const
      {
        return !subset[i];
      }
    };

    struct vertices_not_contained_in
    {
      const vertex_vector& vertices;
      const boost::dynamic_bitset<>& subset;

      vertices_not_contained_in(const vertex_vector& vertices_, const boost::dynamic_bitset<>& subset_)
        : vertices(vertices_),
          subset(subset_)
      {}

      bool operator()(const vertex& v) const
      {
        std::size_t i = &v - &static_cast<const vertex&>(vertices.front());
        return !subset[i];
      }
    };

  public:
    structure_graph() = default;

    structure_graph(vertex_vector vertices, index_type initial_vertex, boost::dynamic_bitset<> exclude)
      : m_vertices(std::move(vertices)),
        m_initial_vertex(initial_vertex),
        m_exclude(std::move(exclude))
    {}

    index_type initial_vertex() const
    {
      return m_initial_vertex;
    }

    std::size_t extent() const
    {
      return m_vertices.size();
    }

    decoration_type decoration(index_type u) const
    {
      return find_vertex(u).decoration;
    }

    std::size_t rank(index_type u) const
    {
      return find_vertex(u).rank;
    }

    const vertex_vector& all_vertices() const
    {
      return m_vertices;
    }

    const std::vector<index_type>& all_predecessors(index_type u) const
    {
      return find_vertex(u).predecessors;
    }

    const std::vector<index_type>& all_successors(index_type u) const
    {
      return find_vertex(u).successors;
    }

    boost::filtered_range<vertices_not_contained_in, const vertex_vector> vertices() const
    {
      return all_vertices() | boost::adaptors::filtered(vertices_not_contained_in(m_vertices, m_exclude));
    }

    boost::filtered_range<integers_not_contained_in, const std::vector<index_type>> predecessors(index_type u) const
    {
      return all_predecessors(u) | boost::adaptors::filtered(integers_not_contained_in(m_exclude));
    }

    boost::filtered_range<integers_not_contained_in, const std::vector<index_type>> successors(index_type u) const
    {
      return all_successors(u) | boost::adaptors::filtered(integers_not_contained_in(m_exclude));
    }

    index_type strategy(index_type u) const
    {
      return find_vertex(u).strategy;
    }

    vertex& find_vertex(index_type u)
    {
      return m_vertices[u];
    }

    const vertex& find_vertex(index_type u) const
    {
      return m_vertices[u];
    }

    const boost::dynamic_bitset<>& exclude() const
    {
      return m_exclude;
    }

    boost::dynamic_bitset<>& exclude()
    {
      return m_exclude;
    }

    bool contains(index_type u) const
    {
      return !m_exclude[u];
    }

    // TODO: avoid this linear time check
    bool is_empty() const
    {
      return m_exclude.all();
    }

    // Returns true if all vertices have a rank and a decoration
    bool is_defined() const
    {
      return std::all_of(m_vertices.begin(), m_vertices.end(), [](const vertex& u) { return u.is_defined(); });
    }

    /// Returns the number of edges in the structure graph.
    std::size_t num_of_edges() {
      std::size_t count = 0;
      for (std::size_t i = 0; i < m_vertices.size(); i++) {
        if (!m_exclude[i]) {
          count += all_successors(i).size();
        }
      }
      return count;
    }
};

template <typename StructureGraph>
std::vector<typename StructureGraph::index_type> structure_graph_predecessors(const StructureGraph& G, typename StructureGraph::index_type u)
{
  std::vector<typename StructureGraph::index_type> result;
  for (auto v: G.predecessors(u))
  {
    result.push_back(v);
  }
  return result;
}

template <typename StructureGraph>
std::vector<typename StructureGraph::index_type> structure_graph_successors(const StructureGraph& G, typename StructureGraph::index_type u)
{
  std::vector<typename StructureGraph::index_type> result;
  for (auto v: G.successors(u))
  {
    result.push_back(v);
  }
  return result;
}

inline
std::ostream& operator<<(std::ostream& out, const structure_graph::decoration_type& decoration)
{
  switch (decoration)
  {
    case structure_graph::d_conjunction : { out << "conjunction"; break; }
    case structure_graph::d_disjunction : { out << "disjunction"; break; }
    case structure_graph::d_true        : { out << "true";        break; }
    case structure_graph::d_false       : { out << "false";       break; }
    default                             : { out << "none";        break; }
  }
  return out;
}

inline
std::ostream& operator<<(std::ostream& out, const structure_graph::vertex& u)
{
  out << "vertex(formula = " << u.formula()
      << ", decoration = " << u.decoration
      << ", rank = " << (u.rank == data::undefined_index() ? std::string("undefined") : std::to_string(u.rank))
      << ", predecessors = " << core::detail::print_list(u.predecessors)
      << ", successors = " << core::detail::print_list(u.successors)
      << ", strategy = " << (u.strategy == undefined_vertex() ? std::string("undefined") : std::to_string(u.strategy))
      << ")";
  return out;
}

template <typename StructureGraph>
std::ostream& print_structure_graph(std::ostream& out, const StructureGraph& G)
{
  auto N = G.all_vertices().size();
  for (std::size_t i = 0; i < N; i++)
  {
    if (G.contains(i))
    {
      const structure_graph::vertex& u = G.find_vertex(i);
      out << std::setw(4) << i << " "
          << "vertex(formula = " << u.formula()
          << ", decoration = " << u.decoration
          << ", rank = " << (u.rank == data::undefined_index() ? std::string("undefined") : std::to_string(u.rank))
          << ", predecessors = " << core::detail::print_list(structure_graph_predecessors(G, i))
          << ", successors = " << core::detail::print_list(structure_graph_successors(G, i))
          << ", strategy = " << (u.strategy == undefined_vertex() ? std::string("undefined") : std::to_string(u.strategy))
          << ")"
          << std::endl;
    }
  }
  if (G.is_empty())
  {
    out << "  empty" << std::endl;
  }
  return out;
}

inline
std::ostream& operator<<(std::ostream& out, const structure_graph& G)
{
  return print_structure_graph(out, G);
}

} // namespace mcrl2::pbes_system



#endif // MCRL2_PBES_STRUCTURE_GRAPH_H
