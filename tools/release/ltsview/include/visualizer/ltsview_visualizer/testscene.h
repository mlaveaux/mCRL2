#ifndef TESTSCENE_H
#define TESTSCENE_H

#include <QOpenGLShaderProgram>
#include <QOpenGLBuffer>
#include <QOpenGLFunctions_4_3_Core>
#include <QOpenGLVertexArrayObject>

#include "glscene.h"
#include <vector>
#include "cluster.h"

#include "glvistree.h"
#include "glscenegraph.h"
#include "glmeshes.h"

#include <map>
#include <tuple>


namespace Primitives{ struct TriangleMesh; }

namespace Test{
// Used in the scenegraph
struct NodeData
{
    int a;
};

// A mesh is a resolution and a meshtype
using Mesh = GlUtil::Meshes::MeshTypes::TriangleMesh;
using MeshRes = std::pair<int, Mesh>; 

struct SceneData{
    int a;
};

using SGNode = GlUtil::PointerTree::PointerTree<NodeData>;
using TestScene = GLScene<NodeData, SceneData>;



class TScene : public TestScene, protected QOpenGLFunctions_4_3_Core{
public:
    explicit TScene(Cluster* root);
    Cluster *m_clusterRoot;
    VisTree::VisTreeNode *m_vistreeRoot;
    // LTSRenderer m_renderer;

    void initializeScene() override;

    void updateScene() override {};

    void renderScene() override;

    void resizeScene(std::size_t width, std::size_t height) override {};

    void rebuildScene() override {};

private:
    QOpenGLShaderProgram* m_prog;
    QOpenGLVertexArrayObject m_vao;
    QOpenGLBuffer m_vbo;
    QOpenGLBuffer m_ibo;
    QOpenGLBuffer m_ssbo;

    GLuint u_n_colors;
    GLuint u_offset_colors;
};
}
#endif