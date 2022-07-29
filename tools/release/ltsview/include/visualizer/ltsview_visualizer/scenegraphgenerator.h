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
#include "vistree.h"
#include "settings.h"
#include <string>
#include <iostream>
#include "cluster.h"

class Cluster;

namespace
{
  struct NodeData
  {
    int model_data_id;
    int model_matrix_id;
  };

  struct ModelData{
      int model_id;
      QMatrix4x4 model_matrix;
  };

  struct Model
  {
    int n_vertices;
    float* vertices;
    int n_triangles;
    int* triangles;
  };

  std::string qvec4ToString(QVector4D vec){
    return "(" + std::to_string(vec.x()) + ", " 
    + std::to_string(vec.y()) + ", " 
    + std::to_string(vec.z()) + ", " 
    + std::to_string(vec.w()) + ")";
  }
} // namespace


struct ConeConvertFunctor{
  QMatrix4x4 current_matrix = QMatrix4x4();
  
  void matrixScale(float x, float y, float z);

  void matrixTranslate(float x, float y, float z);

  void matrixTranslate(const QVector3D& v);

  void matrixRotate(float angle, float x, float y, float z);

  bool topClosed = true;
  float rot = 0;
  VisTree::VisTreeNode* operator()(VisTree::VisTreeNode* parent, Cluster* cluster);

  VisTree::Shape noChildren(Cluster* cluster, bool topClosed, int rot);
};

// TODO: Implement
struct TubeConvertFunctor{
  VisTree::VisTreeNode* operator()(VisTree::VisTreeNode* parent, Cluster* cluster);
};

struct ClusterChildIterator{
  std::vector<Cluster*>::iterator it;
  std::vector<Cluster*>::iterator it_end;
  
  std::vector<Cluster*>::iterator operator++();
  Cluster* operator*();
};

ClusterChildIterator getClusterChildIterator(Cluster* cluster);

typedef SceneGraph<NodeData, ModelData> SG;

template < typename Functor >
VisTree::VisTreeNode* generateClusterTree(Cluster* root){
    Functor f = Functor();
    auto getIterator = std::function<ClusterChildIterator(Cluster*)>(getClusterChildIterator);
    return VisTree::generateVisTree<Cluster, Functor, ClusterChildIterator>(root, f, getIterator);
}
 


class SceneGraphGenerator
{
  public:
    SceneGraphGenerator(Cluster* root) : m_clusterRoot(root) {};

    enum Mode
    {
        TUBES,
        CONES,
    };
    SG* generate(Mode mode);

  private:
    // This should hold all data needed to invoke shaders i.e. verts/tris for general case
    SG* m_sceneGraph;
    // This is generated for specific case
    VisTree::VisTreeNode* m_visTreeRoot;
    Cluster* m_clusterRoot;

    SG* generateTubes(Cluster* root);

    SG* generateCones(Cluster* root);
};

#endif