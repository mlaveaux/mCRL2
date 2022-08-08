// Author(s): Ruben Vink.
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_OPENGL_GLSCENEGRAPH_H
#define MCRL2_OPENGL_GLSCENEGRAPH_H

#include <vector>
#include <QMatrix4x4>
#include <QOpenGLShaderProgram>

template <typename NodeData>
struct SceneGraphNode{
    NodeData data;
    SceneGraphNode* parent;
    std::vector<SceneGraphNode*> children;
};


template <typename NodeData, typename ModelData>
class SceneGraph{
public:
    SceneGraphNode<NodeData>* root;
    ModelData modelData;
    // updates any changed parts of the scenegraph
    void update();

    // completely builds the scenegraph properties
    void rebuild();
};

#endif