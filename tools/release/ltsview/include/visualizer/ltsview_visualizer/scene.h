#ifndef SCENE_H
#define SCENE_H

#include <QOpenGLSHaderProgram>

#include "glscene.h"
#include <vector>
#include "cluster.h"

#include "glvistree.h"
#include "glscenegraph.h"
#include "glmeshes.h"
// #include "renderer.h"


namespace Primitives{ struct TriangleMesh; }


namespace GlLTSView{
    struct NodeData
    {
        int model_id;
        int model_matrix_id;
    };

    using Mesh = GlUtil::MeshTypes::TriangleMesh; 

    struct SceneData{
      std::vector<Mesh> meshes = std::vector<Mesh>();
      std::vector<QMatrix4x4> matrices = std::vector<QMatrix4x4>();
    };

    using SGNode = GlUtil::VectorTree::VectorTree<NodeData>;
    using LTSScene = GLScene<NodeData, SceneData>;

    struct SceneGraphFunctor{
        SceneData& sceneData;
        SceneGraphFunctor(SceneData& sceneData) : sceneData(sceneData) {}

        SGNode* operator()(SGNode* parent,
                                    VisTree::VisTreeNode* cluster);
    };


    class Scene : private LTSScene{
    public:
        Scene(QOpenGLWidget& glwidget, SceneGraph<NodeData, SceneData>& scenegraph, Camera& camera, Cluster* root);
        Cluster *m_clusterRoot;
        VisTree::VisTreeNode *m_vistreeRoot;
        // LTSRenderer m_renderer;

        void initialize() override {};

        void update() override {};

        void render(QPainter& painter) override;

        void resize(std::size_t width, std::size_t height) override {};

        void rebuild() override;

        bool built = false;
    private:
        GLuint m_VAO;
        GLuint m_PositionBuffer;
        GLuint m_NormalBuffer;
        QOpenGLShaderProgram program;
    };
}
#endif