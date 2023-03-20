#ifndef SCENE_H
#define SCENE_H

#include <QOpenGLShaderProgram>
#include <QOpenGLBuffer>
#include <QOpenGLFunctions_4_3_Core>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLExtraFunctions>

#include "glscene.h"
#include <vector>
#include "cluster.h"

#include "glvistree.h"
#include "glscenegraph.h"
#include "glmeshes.h"

#include <map>
#include <tuple>

namespace Primitives
{
struct TriangleMesh;
}

namespace GlLTSView
{
// Used in the shaders
struct VertexData
{
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

struct SceneData
{
  std::vector<NodeData> nodes = std::vector<NodeData>();
  std::vector<Mesh> meshes = std::vector<Mesh>();
  std::vector<GlUtil::Meshes::AABB> aabbs = std::vector<GlUtil::Meshes::AABB>();
  std::vector<QColor> colors = std::vector<QColor>();
  std::vector<QMatrix4x4> matrices = std::vector<QMatrix4x4>();
};

using SGNode = GlUtil::PointerTree::PointerTree<NodeData>;
using LTSScene = GLScene<NodeData, SceneData>;

class ConeDB
{
  public:
  int getCone(GlUtil::Shapes::TruncatedCone* cone, int resolution,
              SceneData& sceneData);

  private:
  typedef std::tuple<qreal, qreal, qreal, int> key_t;
  std::map<key_t, int> cones = std::map<key_t, int>();
};

class SphereDB
{

  public:
  int getSphere(GlUtil::Shapes::Sphere* sphere, int resolution,
                SceneData& sceneData);

  private:
  std::map<std::pair<qreal, int>, int> spheres =
      std::map<std::pair<qreal, int>, int>();
};

template <int sphereRes, int coneRes> struct SceneGraphFunctor
{
  SceneData& sceneData;
  SphereDB sphereDB = SphereDB();
  ConeDB coneDB = ConeDB();

  SceneGraphFunctor(SceneData& sceneData) : sceneData(sceneData)
  {
  }

  SGNode* operator()(SGNode* parent, VisTree::VisTreeNode* cluster);
};

struct Instance
{
  int vertex_offset;
  int vertex_count;
  int index_offset;
  int index_count;
  int instance_count;

  ~Instance()
  {
    vao.destroy();
    vbo.destroy();
    normal_bo.destroy();
    ibo.destroy();
    color_bo.destroy();
    matrix_bo.destroy();
  }

  // VAO
  QOpenGLVertexArrayObject vao;

  // model data
  QOpenGLBuffer vbo;
  QOpenGLBuffer normal_bo;
  QOpenGLBuffer ibo;

  void createModel(Mesh& mesh, int vertex_loc, int normal_loc)
  {
    QOpenGLExtraFunctions* f =
        QOpenGLContext::currentContext()->extraFunctions();
    // VAO
    vao.create();
    vao.bind();

    // IBO
    ibo = QOpenGLBuffer(QOpenGLBuffer::IndexBuffer);
    ibo.create();
    ibo.bind();
    ibo.setUsagePattern(QOpenGLBuffer::UsagePattern::StaticDraw);
    GLuint* indices = new GLuint[3 * mesh.n_triangles];
    int indexindex = 0; // ...
    for (int i = 0; i < mesh.n_triangles; i++)
      for (int j = 0; j < 3; j++)
        indices[indexindex++] = mesh.triangles[i][j];
    ibo.allocate(indices, indexindex * sizeof(GLuint));

    // VBOs
    vbo.create();
    vbo.bind();
    vbo.setUsagePattern(QOpenGLBuffer::UsagePattern::StaticDraw);
    vbo.allocate(&mesh.vertices[0], mesh.vertices.size() * sizeof(GLfloat) * 3);
    f->glVertexAttribPointer(vertex_loc, 3, GL_FLOAT, false, 0, 0);

    normal_bo.create();
    normal_bo.bind();
    normal_bo.setUsagePattern(QOpenGLBuffer::UsagePattern::StaticDraw);
    normal_bo.allocate(&mesh.vertex_normals[0],
                       mesh.vertices.size() * sizeof(GLfloat) * 3);
    f->glVertexAttribPointer(normal_loc, 3, GL_FLOAT, false, 0, 0);
    vao.release();
  }

  // per instance
  std::vector<QMatrix4x4> matrix_cpu;
  std::vector<QVector4D> color_cpu;
  QOpenGLBuffer matrix_bo;
  QOpenGLBuffer color_bo;

  void createInstances(int color_loc, int matrix_loc)
  {
    QOpenGLExtraFunctions* f =
        QOpenGLContext::currentContext()->extraFunctions();
    vao.bind();

    // colors first
    color_bo.create();
    color_bo.bind();
    color_bo.setUsagePattern(QOpenGLBuffer::UsagePattern::DynamicDraw);
    color_bo.allocate(&color_cpu[0], color_cpu.size() * sizeof(GLfloat) * 4);
    f->glVertexAttribPointer(color_loc, 4, GL_FLOAT, false, 0, 0);
    f->glVertexAttribDivisor(color_loc, 1);

    // matrices second
    matrix_bo.create();
    matrix_bo.bind();
    matrix_bo.setUsagePattern(QOpenGLBuffer::UsagePattern::DynamicDraw);
    matrix_bo.allocate(&matrix_cpu[0], matrix_cpu.size() * sizeof(GLfloat) * 4);
    for (int i = 0; i < 4; i++)
    {
      f->glVertexAttribPointer(matrix_loc + i, 4, GL_FLOAT, false,
                               16 * sizeof(GLfloat),
                               (void*)(i * 4 * sizeof(GLfloat)));
      f->glVertexAttribDivisor(matrix_loc + i, 1);
    }
    vao.release();
  }

  void render()
  {
    QOpenGLExtraFunctions* f =
        QOpenGLContext::currentContext()->extraFunctions();
    vao.bind();
    f->glDrawElementsInstanced(GL_TRIANGLES, index_count, GL_UNSIGNED_INT, 0,
                               instance_count);
    vao.release();
  }
};

class Scene : public LTSScene, protected QOpenGLFunctions_4_3_Core
{
  public:
  explicit Scene(QOpenGLWidget* glwidget, Cluster* root);
  Cluster* m_clusterRoot;
  VisTree::VisTreeNode* m_vistreeRoot;
  // LTSRenderer m_renderer;

  void initializeScene() override;

  void updateScene() override{};

  void renderScene() override;

  void resizeScene(std::size_t width, std::size_t height) override{};

  void rebuildScene() override;

  private:
  bool built = false;
  bool buffers_exist = false;
  int m_ssbo_count = 0;

  std::vector<Instance*> m_instances;

  GLuint m_uniform_view_loc;
  GLuint m_uniform_proj_loc;
  GLuint m_uniform_alpha_loc;

  // helper functions
  bool reportOpenGLError();

  QOpenGLShaderProgram m_prog;

  std::vector<VertexData> m_vertexData;
  std::vector<QVector3D> m_vertices;
  std::vector<QVector3D> m_normals;
  std::vector<QColor> m_colors;
  std::vector<GlUtil::Triangle> m_triangles;
};
} // namespace GlLTSView
#endif