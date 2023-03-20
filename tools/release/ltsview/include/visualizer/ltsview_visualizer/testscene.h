#ifndef TESTSCENE_H
#define TESTSCENE_H

#include <QOpenGLShaderProgram>
#include <QOpenGLBuffer>
#include <QOpenGLExtraFunctions>
#include <QOpenGLFunctions_3_3_Core>
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


struct Instance
{
  int vertex_offset;
  int vertex_count;
  int index_offset;
  int index_count;
  int instance_count;

  ~Instance()
  {
    vbo.destroy();
    normal_bo.destroy();
    ibo.destroy();
    color_bo.destroy();
    matrix_bo.destroy();
  }


  // model data
  QOpenGLBuffer vbo;
  QOpenGLBuffer normal_bo;
  QOpenGLBuffer ibo;

  void createModel(QOpenGLFunctions_3_3_Core* f, Mesh& mesh);

  // per instance
  std::vector<QMatrix4x4> matrix_cpu;
  std::vector<QVector4D> color_cpu;
  QOpenGLBuffer matrix_bo;
  QOpenGLBuffer color_bo;

  // locations
  int m_vertex_loc;
  int m_normal_loc;
  int m_color_loc;
  int m_matrix_loc;

  void createInstances(QOpenGLFunctions_3_3_Core* f);

  void render(QOpenGLFunctions_3_3_Core* f, QOpenGLShaderProgram& shader);
};


class TScene : public TestScene, private QOpenGLFunctions_3_3_Core{
public:
    explicit TScene(Cluster* root);
    Cluster *m_clusterRoot;
    VisTree::VisTreeNode *m_vistreeRoot;
    // LTSRenderer m_renderer;

    void initializeScene() override;

    void updateScene() override {};

    void renderScene() override;

    void resizeScene(std::size_t width, std::size_t height) override
    {
      m_camera->setViewport(width, height);
    };

    void rebuildScene() override;


private:
    std::vector<Instance*> m_instances;
    QOpenGLVertexArrayObject m_vao;
    QOpenGLShaderProgram m_prog;


    int m_total_vert = -1;
    int m_total_tris = -1;
    
    bool m_cpu_allocated = false;

    GLuint u_n_colors;
    GLuint u_offset_colors;

    int m_vertex_loc;
    int m_normal_loc;
    int m_color_loc;
    int m_model_matrix_loc;

    bool m_built = false;
};
}
#endif