// Author(s): Bas Ploeger
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <list>
#include "rtree.h"

class Rectangle
{
  public:
    Rectangle() {}

    Rectangle(const QVector2D& lowc, const QVector2D& highc):
      low_corner(lowc), high_corner(highc)
    { }

    bool contains(const QVector2D& point) const;

    void ensureMinSideLength(float min_side_length);

    float minDistanceTo(const QVector2D& point) const;


    // We assume the following invariant for any rectangle:
    // low_corner.x() <= high_corner.x() && low_corner.y() <= high_corner.y()
    QVector2D low_corner;
    QVector2D high_corner;
};

bool Rectangle::contains(const QVector2D& point) const
{
  return point.x() >= low_corner.x() && high_corner.x() >= point.x() &&
         point.y() >= low_corner.y() && high_corner.y() >= point.y();
}

void Rectangle::ensureMinSideLength(float min_side_length)
{
  float lc_x = low_corner.x();
  float lc_y = low_corner.y();
  float hc_x = high_corner.x();
  float hc_y = high_corner.y();
  float diff_x = hc_x - lc_x;
  if (diff_x < min_side_length)
  {
    hc_x += 0.5f * (min_side_length - diff_x);
    lc_x -= 0.5f * (min_side_length - diff_x);
  }
  float diff_y = hc_y - lc_y;
  if (diff_y < min_side_length)
  {
    hc_y += 0.5f * (min_side_length - diff_y);
    lc_y -= 0.5f * (min_side_length - diff_y);
  }
  low_corner = QVector2D(lc_x, lc_y);
  high_corner = QVector2D(hc_x, hc_y);
}


class RNode
{
  public:
    virtual ~RNode()
    { }

    Rectangle boundingBox() const
    {
      return bounding_box;
    }

    float centerX() const
    {
      return bounding_box.low_corner.x() + 0.5f * (bounding_box.high_corner.x()
             - bounding_box.low_corner.x());
    }

    float centerY() const
    {
      return bounding_box.low_corner.y() + 0.5f * (bounding_box.high_corner.y()
             - bounding_box.low_corner.y());
    }

    virtual void computeBoundingBox()
    { }

    virtual void deletePoint(const QVector2D&)
    { }

    virtual bool hasChildren() const
    {
      return false;
    }

  protected:
    Rectangle bounding_box;

    static const float MIN_BB_SIDE_LENGTH;
};

const float RNode::MIN_BB_SIDE_LENGTH = 0.001f;

static inline
bool lessThanCenterX(RNode* a, RNode* b)
{
  // See cluster.cpp for an explanation on the a != b clause.
  return a != b && a->centerX() < b->centerX();
}

static inline
bool lessThanCenterY(RNode* a, RNode* b)
{
  // See cluster.cpp for an explanation on the a != b clause.
  return a != b && a->centerY() < b->centerY();
}


class RTreeNode: public RNode
{
  public:
    RTreeNode()
    { }

    ~RTreeNode()
    {
      for (std::list< RNode* >::iterator ci = children.begin();
           ci != children.end(); ++ci)
      {
        delete *ci;
      }
    }

    void computeBoundingBox();

    void deletePoint(const QVector2D& point);

    bool hasChildren() const
    {
      return !children.empty();
    }

    std::list< RNode* > children;
};

void RTreeNode::computeBoundingBox()
{
  float lc_x = std::numeric_limits< float >::max();
  float lc_y = std::numeric_limits< float >::max();
  float hc_x = std::numeric_limits< float >::min();
  float hc_y = std::numeric_limits< float >::min();
  std::list< RNode* >::iterator ci;
  for (ci = children.begin(); ci != children.end(); ++ci)
  {
    Rectangle child_bb = (**ci).boundingBox();
    lc_x = std::min(lc_x, child_bb.low_corner.x());
    lc_y = std::min(lc_y, child_bb.low_corner.y());
    hc_x = std::max(hc_x, child_bb.high_corner.x());
    hc_y = std::max(hc_y, child_bb.high_corner.y());
  }
  bounding_box = Rectangle(QVector2D(lc_x, lc_y), QVector2D(hc_x, hc_y));
  bounding_box.ensureMinSideLength(MIN_BB_SIDE_LENGTH);
}

void RTreeNode::deletePoint(const QVector2D& point)
{
  std::list< RNode* >::iterator ci = children.begin();
  while (ci != children.end())
  {
    RNode* child = *ci;
    if (child->boundingBox().contains(point))
    {
      child->deletePoint(point);
      if (child->hasChildren())
      {
        child->computeBoundingBox();
      }
      else
      {
        delete child;
        ci = children.erase(ci);
        continue;
      }
    }
    ++ci;
  }
}


class RTreeLeaf: public RNode
{
  public:
    RTreeLeaf(const QVector2D& p):
      point(p)
    { }

    ~RTreeLeaf()
    { }

    void computeBoundingBox();

    bool hasChildren() const
    {
      return false;
    }

    void deletePoint(const QVector2D&)
    { }

    QVector2D point;
};

void RTreeLeaf::computeBoundingBox()
{
  float lc_x = point.x();
  float lc_y = point.y();
  float hc_x = point.x();
  float hc_y = point.y();
  bounding_box = Rectangle(QVector2D(lc_x, lc_y), QVector2D(hc_x, hc_y));
  bounding_box.ensureMinSideLength(MIN_BB_SIDE_LENGTH);
}

struct QueueElement
{
  RNode* node;
  float distance;
  QueueElement(RNode* n, float d): node(n), distance(d) { }
};

class QueueElementLessThan
{
  public:
    bool operator()(const QueueElement& q, const QueueElement& r)
    {
      return q.distance < r.distance;
    }
};

class QueueElementGreaterThan
{
  public:
    bool operator()(const QueueElement& q, const QueueElement& r)
    {
      return q.distance > r.distance;
    }
};

class MinRectPointDistance
{
  public:
    float operator()(const Rectangle& rectangle, const QVector2D& point)
    {
      if (rectangle.contains(point))
      {
        return 0.0f;
      }
      QVector2D diff_low = rectangle.low_corner - point;
      QVector2D diff_high = rectangle.high_corner - point;
      return std::min(diff_low.x() * diff_low.x(), diff_high.x() * diff_high.x()) +
             std::min(diff_low.y() * diff_low.y(), diff_high.y() * diff_high.y());
    }
};

class MaxRectPointDistance
{
  public:
    float operator()(const Rectangle& rectangle, const QVector2D& point)
    {
      QVector2D diff_low = rectangle.low_corner - point;
      QVector2D diff_high = rectangle.high_corner - point;
      return std::max(diff_low.x() * diff_low.x(), diff_high.x() * diff_high.x()) +
             std::max(diff_low.y() * diff_low.y(), diff_high.y() * diff_high.y());
    }
};

template< class Distance, class Compare >
class NeighbourFinder
{
  public:
    typedef std::priority_queue< QueueElement, std::vector< QueueElement >,
            Compare > PriorityQueue;

    void startNewSearch(RNode* root, const QVector2D& point)
    {
      queue = PriorityQueue();
      has_found_neighbour = false;
      query_point = point;
      queue.push(QueueElement(root, 0.0f));
    }

    void findNextNeighbour()
    {
      has_found_neighbour = false;
      while (!queue.empty())
      {
        RNode* top_node = queue.top().node;
        queue.pop();
        RTreeLeaf* leaf = dynamic_cast< RTreeLeaf* >(top_node);
        if (leaf != NULL)
        {
          has_found_neighbour = true;
          found_neighbour = leaf->point;
          return;
        }
        RTreeNode* node = dynamic_cast< RTreeNode* >(top_node);
        std::list< RNode* >::iterator ci;
        for (ci = node->children.begin(); ci != node->children.end(); ++ci)
        {
          RNode* child = *ci;
          queue.push(QueueElement(child, Distance()(child->boundingBox(),
                                  query_point)));
        }
      }
    }

    bool hasFoundNeighbour() const
    {
      return has_found_neighbour;
    }

    QVector2D foundNeighbour() const
    {
      return found_neighbour;
    }

  private:
    QVector2D found_neighbour;
    bool has_found_neighbour;
    QVector2D query_point;
    PriorityQueue queue;
};

RTree::~RTree()
{
  if (root != NULL)
  {
    delete root;
  }
}

void RTree::deletePoint(const QVector2D& point)
{
  root->deletePoint(point);
  if (!root->hasChildren())
  {
    delete root;
    root = NULL;
  }
}

void RTree::findFarthestNeighbour(const QVector2D& point)
{
  NeighbourFinder< MaxRectPointDistance, QueueElementLessThan > finder;
  finder.startNewSearch(root, point);
  finder.findNextNeighbour();
  neighbour_found = finder.hasFoundNeighbour();
  if (neighbour_found)
  {
    neighbour = finder.foundNeighbour();
  }
}

void RTree::findNearestNeighbour(const QVector2D& point)
{
  NeighbourFinder< MinRectPointDistance, QueueElementGreaterThan > finder;
  finder.startNewSearch(root, point);
  finder.findNextNeighbour();
  neighbour_found = finder.hasFoundNeighbour();
  if (neighbour_found)
  {
    neighbour = finder.foundNeighbour();
  }
}

// The maximum number of children for a node in the RTree
const unsigned int PackedRTreeBuilder::MAX_FANOUT = 50;

void PackedRTreeBuilder::buildRTree()
{
  // This algorithm for constructing a packed R-tree from a set of spatial
  // objects (in our case just points), is called Sort-Tile-Recursive and is
  // described in:
  // "STR: a simple and efficient algorithm for R-tree packing", S.T.
  // Leutenegger, M.A. Lopez and J. Edgington, Technical Report TR-97-14,
  // Institute for Computer Applications in Science and Engineering, 1997.
  std::vector< RNode* > roots;
  roots.reserve(points.size());
  std::vector< QVector2D >::iterator pi;
  for (pi = points.begin(); pi != points.end(); ++pi)
  {
    RTreeLeaf* leaf = new RTreeLeaf(*pi);
    leaf->computeBoundingBox();
    roots.push_back(leaf);
  }
  while (roots.size() > 1)
  {
    std::vector< RNode* > new_roots;
    unsigned int slice_size = static_cast<unsigned int>(MAX_FANOUT * std::ceil(std::sqrt(std::ceil(
                                static_cast< double >(roots.size()) /
                                static_cast< double >(MAX_FANOUT)))));
    sort(roots.begin(), roots.end(), lessThanCenterX);
    for (unsigned int i = 0; i < roots.size(); i += slice_size)
    {
      unsigned int this_slice_size = std::min(slice_size,
                                              static_cast<unsigned int>(roots.size()) - i);
      std::vector< RNode* >::iterator slice_begin = roots.begin() + i;
      std::vector< RNode* >::iterator slice_end = slice_begin + this_slice_size;
      sort(slice_begin, slice_end, lessThanCenterY);
      for (unsigned int j = 0; j < this_slice_size; j += MAX_FANOUT)
      {
        unsigned int num_children = std::min(MAX_FANOUT, this_slice_size - j);
        std::vector< RNode* >::iterator children_begin = slice_begin + j;
        RTreeNode* node = new RTreeNode();
        node->children.assign(children_begin, children_begin + num_children);
        node->computeBoundingBox();
        new_roots.push_back(node);
      }
    }
    roots.swap(new_roots);
  }
  rtree = new RTree(roots[0]);
}
