// Author(s): Ruben Vink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef MCRL2_OPENGL_SCENEGRAPHGENERATOR_H
#define MCRL2_OPENGL_SCENEGRAPHGENERATOR_H

#include <vector>
#include <QVector3D>
#include <QMatrix4x4>
#include <stack>
#include "glscenegraph.h"
#include "glvistree.h"
#include "settings.h"
#include <string>
#include <iostream>
#include "cluster.h"
#include <QtCore>

class Cluster;

namespace
{
std::string qvec4ToString(QVector4D vec)
{
  return "(" + std::to_string(vec.x()) + ", " + std::to_string(vec.y()) + ", " +
         std::to_string(vec.z()) + ", " + std::to_string(vec.w()) + ")";
}
} // namespace

struct ConeConvertFunctor
{
  QMatrix4x4 current_matrix = QMatrix4x4();

  void matrixScale(float x, float y, float z);

  void matrixTranslate(float x, float y, float z);

  void matrixTranslate(const QVector3D& v);

  void matrixRotate(float angle, float x, float y, float z);

  bool topClosed = true;
  qreal max_rank = -1;
  QColor c_top, c_bot;
  VisTree::VisTreeNode* operator()(VisTree::VisTreeNode* parent,
                                   Cluster* cluster);
};

// TODO: Implement
struct TubeConvertFunctor
{
  qreal max_rank = -1;
  VisTree::VisTreeNode* operator()(VisTree::VisTreeNode* parent,
                                   Cluster* cluster);
};

template <typename Functor>
VisTree::VisTreeNode* generateClusterTree(Cluster* root, Functor& f)
{
  std::function<std::vector<Cluster*>::iterator(Cluster*)> getChildBegin =
      [](Cluster* node) { return node->descendants.begin(); };
  std::function<std::vector<Cluster*>::iterator(Cluster*)> getChildEnd =
      [](Cluster* node) { return node->descendants.end(); };
  return VisTree::generateVisTree<Cluster, Functor,
                                  std::vector<Cluster*>::iterator>(
      root, f, getChildBegin, getChildEnd);
}

class VisTreeGenerator
{
  public:
  enum Mode
  {
    TUBES,
    CONES,
  };
  static VisTree::VisTreeNode* generate(Mode mode, Cluster* root);

  private:
  // static class so private constructor
  VisTreeGenerator()
  {
  }

  static VisTree::VisTreeNode* generateTubes(Cluster* root, int max_rank);

  static VisTree::VisTreeNode* generateCones(Cluster* root, int max_rank);
};

#endif