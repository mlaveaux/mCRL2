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

namespace{
    template <typename NodeData>
    struct SceneGraphNode{
        NodeData data;
        SceneGraphNode* parent;
        std::vector<SceneGraphNode*> children;
    };
}

struct AbstractSceneModelData{
    uint dirty = 0; // Use bit flags to set what part is dirty. Implementation specific
    virtual void toGPU() = 0;
};

template < typename T >
struct SceneModelData : AbstractSceneModelData{
    void toGPU(QOpenGLShaderProgram& shader_program){
        if (dirty != 0){
            fix();
            dirty = 0;
        }
    }
    void fix();
    T data;
};

template <typename NodeData, typename ModelData>
class SceneGraph{
public:
    SceneGraphNode<NodeData>* root;
    SceneModelData<ModelData> modelData;
    // updates any changed parts of the scenegraph
    void update();

    // completely builds the scenegraph properties
    void rebuild();
};

#endif