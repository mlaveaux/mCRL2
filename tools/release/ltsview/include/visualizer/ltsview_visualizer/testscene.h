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
// Used in the shaders
struct VertexData{
    int model_index;
    int vertex_index;
};

// Used in the scenegraph
struct NodeData
{
    int model_id;
    int data_id;
};

// A mesh is a resolution and a meshtype
using Mesh = GlUtil::Meshes::MeshTypes::TriangleMesh;
using MeshRes = std::pair<int, Mesh>; 

struct SceneData{
    std::vector<NodeData> nodes = std::vector<NodeData>();
    std::vector<Mesh> meshes = std::vector<Mesh>();
    std::vector<GlUtil::Meshes::AABB> aabbs = std::vector<GlUtil::Meshes::AABB>();
    std::vector<QColor> colors = std::vector<QColor>();
    std::vector<QMatrix4x4> matrices = std::vector<QMatrix4x4>();
    void clear();
};

using SGNode = GlUtil::PointerTree::PointerTree<NodeData>;
using TestScene = GLScene<NodeData, SceneData>;

// Scenegraph generator classes
class ConeDB{
    public:
        int getCone(GlUtil::Shapes::TruncatedCone* cone, int resolution, SceneData& sceneData);

    private:
        typedef std::tuple<qreal, qreal, qreal, int> key_t;
        std::map<key_t, int> cones = std::map<key_t, int>(); 
};

class SphereDB{
    
    public:
        int getSphere(GlUtil::Shapes::Sphere* sphere, int resolution, SceneData& sceneData);

    private:
        std::map<std::pair<qreal, int>, int> spheres = std::map<std::pair<qreal, int>, int>();
};

template < int sphereRes, int coneRes >
struct SceneGraphFunctor{
    SceneData& sceneData;
    SphereDB sphereDB = SphereDB();
    ConeDB coneDB = ConeDB();
    
    SceneGraphFunctor(SceneData& sceneData) : sceneData(sceneData) {}

    SGNode* operator()(SGNode* parent,
                                VisTree::VisTreeNode* cluster);
};
// End scenegraph generator classes



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

    void rebuildScene() override;


private:
    int m_ssbo_count = 0;
    void createBufferObject(QOpenGLBuffer& buff, const char* buff_name, const void* data, std::size_t size, QOpenGLExtraFunctions *f, QOpenGLFunctions_4_3_Core *f430);
    
    QOpenGLShaderProgram m_prog;
    QOpenGLVertexArrayObject m_vao;
    QOpenGLBuffer m_vbo;
    QOpenGLBuffer m_ibo;
    QOpenGLBuffer m_ssbo1;
    QOpenGLBuffer m_ssbo2;
    QOpenGLBuffer m_ssbo3;
    QOpenGLBuffer m_ssbo4;

    int m_total_vert = -1;
    int m_total_tris = -1;

    int* m_vertexdata_cpu;
    float* m_vertices_cpu;
    float* m_normals_cpu;
    int* m_triangles_cpu;
    float* m_matrices_cpu;
    float* m_colors_cpu;
    
    bool m_cpu_allocated = false;
    void reallocateStorage();

    GLuint u_n_colors;
    GLuint u_offset_colors;

    bool m_built = false;
};
}
#endif