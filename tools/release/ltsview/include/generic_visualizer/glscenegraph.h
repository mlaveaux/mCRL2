// Author(s): Ruben Vink.
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_OPENGL_SCENEGRAPH_H
#define MCRL2_OPENGL_SCENEGRAPH_H

#include <vector>
#include <QMatrix4x4>
#include <QOpenGLShaderProgram>

#include "glvectortree.h"



template <typename NodeData, typename SceneData> class SceneGraph
{
  public:
  GlUtil::PointerTree::PointerTree<NodeData>* root;
  SceneData sceneData;
  // updates any changed parts of the scenegraph
  void update();

  // completely builds the scenegraph properties
  void rebuild();
};

#endif