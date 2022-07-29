// Author(s): Bas Ploeger
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef RTREE_H
#define RTREE_H

#include <vector>
#include "vectors.h"

class RNode;

class RTree
{
  public:
    RTree():
      root(NULL),
      neighbour_found(false)
    { }

    RTree(RNode* root_):
      root(root_),
      neighbour_found(false)
    { }

    ~RTree();
    void deletePoint(const QVector2D& point);
    void findNearestNeighbour(const QVector2D& point);
    void findFarthestNeighbour(const QVector2D& point);

    bool hasFoundNeighbour() const
    {
      return neighbour_found;
    }

    QVector2D foundNeighbour() const
    {
      return neighbour;
    }

  private:
    RNode* root;
    bool neighbour_found;
    QVector2D neighbour;
};


class PackedRTreeBuilder
{
  public:
    PackedRTreeBuilder():
      rtree(NULL)
    { }

    ~PackedRTreeBuilder()
    { }

    void addPoint(const QVector2D& point)
    {
      points.push_back(point);
    }

    void buildRTree();

    RTree* getRTree() const
    {
      return rtree;
    }

  private:
    std::vector< QVector2D > points;
    RTree* rtree;
    // The maximum number of children for a node in the RTree
    static const unsigned int MAX_FANOUT;
};

#endif
