#include "testscene.h"
#include "cluster.h"
#include "vistreegenerator.h"
#include "settings.h"
// Generic opengl files
#include "glutil.h"
#include "glscenegraph.h"
#include "glprimitivefactories.h"
#include "arcballcamera.h"

// std includes
#include <fstream>
#include <sstream>
#include <chrono>

#include <QOpenGLFunctions>
#include <QOpenGLExtraFunctions>

#undef DEBUG_PRINT_VP
#undef DEBUG_PRINT_MATRICES
#undef DEBUG_PRINT_VERTICES

inline int cone_encode_resolution(int res, bool a, bool b)
{
  return res | (a ? (1 << 30) : 0) | (b ? (1 << 31) : 0);
}
int Test::ConeDB::getCone(GlUtil::Shapes::TruncatedCone* cone, int resolution,
                          SceneData& sceneData)
{
  key_t key = {
      cone->radius_top, cone->radius_bot, cone->height,
      cone_encode_resolution(resolution, cone->fill_top, cone->fill_bot)};
  auto it = cones.find(key);
  if (it != cones.end())
  {
    return it->second;
  }
  else
  {
    // create the mesh
    Test::Mesh mesh =
        GlUtil::DefaultFactories::TruncatedConeFactory<Mesh>::createPrimitive(
            cone, resolution);
    int mesh_id = (int)sceneData.meshes.size();
    sceneData.meshes.push_back(mesh);
    cones[key] = mesh_id;
    return mesh_id;
  }
}

int Test::SphereDB::getSphere(GlUtil::Shapes::Sphere* sphere, int resolution,
                              SceneData& sceneData)
{
  auto it = spheres.find({sphere->radius, resolution});
  if (it != spheres.end())
  {
    return it->second;
  }
  else
  {
    // create the mesh
    Test::Mesh mesh =
        GlUtil::DefaultFactories::SphereFactory<Mesh>::createPrimitive(
            sphere, resolution);
    int mesh_id = (int)sceneData.meshes.size();
    sceneData.meshes.push_back(mesh);
    spheres[{sphere->radius, resolution}] = mesh_id;
    return mesh_id;
  }
}

void Test::Instance::createModel(QOpenGLFunctions_3_3_Core* f, Mesh& mesh)
{
  mCRL2log(mcrl2::log::debug) << "Instance::createModel(...): " << std::endl;
  std::string yeet = "test.ply";
  GlUtil::Meshes::Export::exportMesh<GlUtil::Meshes::MeshTypes::TriangleMesh, GlUtil::Meshes::Export::FileType::PLY>(yeet, mesh);

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
  float* temp_vertices = new float[mesh.n_vertices * 3];
  int temp_vertex_ind = 0;
  for (auto& vec : mesh.vertices)
  {
    temp_vertices[temp_vertex_ind++] = vec.x();
    temp_vertices[temp_vertex_ind++] = vec.y();
    temp_vertices[temp_vertex_ind++] = vec.z();
  }
  vbo.allocate(temp_vertices, temp_vertex_ind * sizeof(GLfloat));

  normal_bo.create();
  normal_bo.bind();
  normal_bo.setUsagePattern(QOpenGLBuffer::UsagePattern::StaticDraw);
  normal_bo.allocate(&mesh.vertex_normals[0],
                     mesh.vertices.size() * sizeof(GLfloat) * 3);
}

void Test::Instance::createInstances(QOpenGLFunctions_3_3_Core* f)
{
  mCRL2log(mcrl2::log::debug) << "Openglcontext: " << f << " inst  count " << instance_count << std::endl;

  glCheckError();
  // colors first
  color_bo.create();
  color_bo.bind();
  color_bo.setUsagePattern(QOpenGLBuffer::UsagePattern::DynamicDraw);
  float* temp_color = new float[color_cpu.size() * 4];
  int temp_color_ind = 0;
  for (auto& c : color_cpu)
  {
    temp_color[temp_color_ind++] = c.x();
    temp_color[temp_color_ind++] = c.y();
    temp_color[temp_color_ind++] = c.z();
    temp_color[temp_color_ind++] = c.w();
  }
  color_bo.allocate(temp_color, instance_count * 4 * sizeof(GLfloat));
  // matrices second
  matrix_bo.create();
  matrix_bo.bind();
  matrix_bo.setUsagePattern(QOpenGLBuffer::UsagePattern::DynamicDraw);
  float* temp_matrix = new float[matrix_cpu.size() * 16];
  int temp_matrix_ind = 0;
  for (auto& m : matrix_cpu)
  {
    float* temp = new float[16];
    m.transposed().copyDataTo(temp);
    for (int i{0}; i < 16; ++i)
    {
      temp_matrix[temp_matrix_ind++] = temp[i];
    }
  }
  matrix_bo.allocate(temp_matrix, instance_count * 16 * sizeof(GLfloat));
  glCheckError();
}

void Test::Instance::render(QOpenGLFunctions_3_3_Core* f, QOpenGLShaderProgram& shader)
{
  //mCRL2log(mcrl2::log::debug) << "Openglcontext: " << f << std::endl;

  //glCheckError();
  //int bound_vao;
  //f->glGetIntegerv(GL_VERTEX_ARRAY_BINDING, &bound_vao);
  //glCheckError();
  //mCRL2log(mcrl2::log::debug) << "Instance::render bound_vao: " << bound_vao << std::endl;

  m_vertex_loc = f->glGetAttribLocation(shader.programId(), "vertex");
  m_normal_loc = f->glGetAttribLocation(shader.programId(), "normal");
  m_color_loc = f->glGetAttribLocation(shader.programId(), "color");
  m_matrix_loc = f->glGetAttribLocation(shader.programId(), "model_mat");

  vbo.bind();
  glCheckError();
  f->glVertexAttribPointer(m_vertex_loc, 3, GL_FLOAT, false, 0, 0);
  glCheckError();
  normal_bo.bind();
  shader.setAttributeBuffer(m_normal_loc, GL_FLOAT, 0, 3);
  glCheckError();
  color_bo.bind();
  f->glVertexAttribPointer(m_color_loc, 4, GL_FLOAT, false, 0, 0);
  glCheckError();
  matrix_bo.bind();
  f->glVertexAttribPointer(m_matrix_loc, 4, GL_FLOAT, GL_FALSE, sizeof(GLfloat) * 4 * 4,
                        (void*)(0));
  f->glVertexAttribPointer(m_matrix_loc+1, 4, GL_FLOAT, GL_FALSE, sizeof(GLfloat) * 4 * 4,
                        (void*)(sizeof(float) * 4));
  f->glVertexAttribPointer(m_matrix_loc + 2, 4, GL_FLOAT, GL_FALSE,
                           sizeof(GLfloat) * 4 * 4,
                        (void*)(sizeof(float) * 8));
  f->glVertexAttribPointer(m_matrix_loc + 3, 4, GL_FLOAT, GL_FALSE,
                           sizeof(GLfloat) * 4 * 4,
                        (void*)(sizeof(float) * 12));
  glCheckError();

  ibo.bind();
  f->glDrawElementsInstanced(GL_TRIANGLES, 3*index_count, GL_UNSIGNED_INT, 0,
                             instance_count);
  glCheckError();
}

void Test::SceneData::clear()
{
  nodes.resize(0);

  meshes.resize(0);

  aabbs.resize(0);

  colors.resize(0);

  matrices.resize(0);
}

/// TODO: Fix hardcoded resolutions
template <int sphereRes, int coneRes>
Test::SGNode* Test::SceneGraphFunctor<sphereRes, coneRes>::operator()(
    SGNode* parent, VisTree::VisTreeNode* vistreenode)
{
  SGNode* node = new SGNode();
  Mesh mesh;
  node->data.data_id = sceneData.matrices.size();
  sceneData.matrices.push_back(vistreenode->data.matrix);
  sceneData.colors.push_back(vistreenode->data.color);
  switch (vistreenode->data.shape->getShapeType())
  {
  case GlUtil::ShapeType::SPHERE:
    node->data.model_id = sphereDB.getSphere(
        static_cast<GlUtil::Shapes::Sphere*>(vistreenode->data.shape),
        sphereRes, sceneData);
    break;
  case GlUtil::ShapeType::TRUNCATED_CONE:
    node->data.model_id = coneDB.getCone(
        static_cast<GlUtil::Shapes::TruncatedCone*>(vistreenode->data.shape),
        std::max((1 + (coneRes / vistreenode->num_children)) *
                     vistreenode->num_children,
                 6 * vistreenode->num_children) +
            1,
        sceneData);
    break;
  default:
    /// TODO: Log unknown shape
    node->data.model_id = -1;
    node->data.data_id = -1;
    break;
  }
  sceneData.nodes.push_back(node->data);
  return node;
}

const char* VS =
    "#version 330\n"
    "// per vertex\n"
    "layout(location = 0) in vec3 vertex;\n"
    "layout(location = 1) in vec3 normal;\n"

    "// per shape\n"
    "layout(location = 2) in vec4 color;\n"
    "layout(location = 3) in mat4 model_mat;\n"

    "uniform mat4 u_view;\n"
    "uniform mat4 u_proj;\n"
    "uniform float u_alpha;\n"

    "out vec4 vColor;\n"
    "out vec4 vNormal;\n"

    "void main(void)\n"
    "{\n"
    "  vColor = vec4(color.rgb, u_alpha);\n"
    "  vNormal = u_proj * u_view * model_mat * vec4(normal, 0);\n"
    "  gl_Position = u_proj * u_view * model_mat * vec4(vertex, 1);\n"
    "}";

const char* FS = "#version 330\n"
                 "in vec4 vNormal;\n"
                 "in vec4 vColor;\n"
                 "out vec4 fColor;\n"
                 "void main(void)\n"
                 "{\n"
                 "   fColor = vColor;\n"
                 "}";

Test::TScene::TScene(Cluster* root) : TestScene()
{
  SceneGraph<Test::NodeData, Test::SceneData> sg =
      SceneGraph<Test::NodeData, Test::SceneData>();
  m_camera = new ArcballCamera();
  m_camera->reset();
  m_scenegraph = sg;
  m_clusterRoot = root;
}

static QVector4D color_to_vector(const QColor& color)
{
  return {static_cast<float>(color.redF()), static_cast<float>(color.greenF()),
          static_cast<float>(color.blueF()),
          static_cast<float>(color.alphaF())};
}

void Test::TScene::rebuildScene()
{
  m_total_tris = 0;
  m_total_vert = 0;
  m_prog.bind();
  std::cout << "Scene::rebuild() called" << std::endl;
  if (!m_clusterRoot)
  {
    std::cout << "No root found." << std::endl;
    return;
  }
  auto start = std::chrono::high_resolution_clock::now();
  // clear the scene data
  m_scenegraph.sceneData.meshes.clear();
  m_scenegraph.sceneData.nodes.clear();
  m_instances.clear();
  if (m_built)
    free(m_scenegraph.root);
  if (m_built)
    free(m_vistreeRoot);
  if (m_built)
    m_scenegraph.sceneData = SceneData();
  m_vertex_loc = m_prog.attributeLocation("vertex");
  m_normal_loc = m_prog.attributeLocation("normal");
  m_color_loc = m_prog.attributeLocation("color");
  m_model_matrix_loc = m_prog.attributeLocation("model_mat");



  /// TODO: Change settings.h to have enum of modes instead of bool
  mCRL2log(mcrl2::log::debug)
      << "----------------------------------------------------------"
      << (Settings::instance().clusterVisStyleTubes.value() ? "TUBES" : "CONES")
      << std::endl;
  m_vistreeRoot = VisTreeGenerator::generate(
      Settings::instance().clusterVisStyleTubes.value()
          ? VisTreeGenerator::Mode::TUBES
          : VisTreeGenerator::Mode::CONES,
      m_clusterRoot);
  const int sphereRes = 3;
  const int coneRes = 50;
  SceneGraphFunctor<sphereRes, coneRes> sgf(m_scenegraph.sceneData);

  m_scenegraph.root =
      GlUtil::fold_tree<SGNode, VisTree::VisTreeNode,
                        Test::SceneGraphFunctor<sphereRes, coneRes>,
                        VisTree::VisTreeNode::childIterator>(
          nullptr, m_vistreeRoot, SGNode::setParent, SGNode::addChild, sgf,
          VisTree::VisTreeNode::getChildBegin,
          VisTree::VisTreeNode::getChildEnd);

  int n_verts = 0;
  int n_tris = 0;
  for (auto& mesh : m_scenegraph.sceneData.meshes)
  {
    m_instances.emplace_back(new Instance());
    Instance& new_instance = *m_instances.back();
    new_instance.vertex_count = mesh.n_vertices;
    new_instance.index_count = mesh.n_triangles;
    new_instance.createModel(this, mesh);

    n_verts += mesh.n_vertices;
    n_tris += mesh.n_triangles;
  }
  for (auto& node : m_scenegraph.sceneData.nodes)
  {
    Instance& instance = *m_instances[node.model_id];
    instance.color_cpu.push_back(
        color_to_vector(m_scenegraph.sceneData.colors[node.data_id]));
    instance.matrix_cpu.push_back(
        m_scenegraph.sceneData.matrices[node.data_id]);
    instance.instance_count++;
  }

  for (auto& instance : m_instances)
  {
    instance->createInstances(this);
    m_total_vert += instance->vertex_count * instance->instance_count;
    m_total_tris += instance->index_count * instance->instance_count;
  }

  auto end = std::chrono::high_resolution_clock::now();

  mCRL2log(mcrl2::log::debug)
      << "Took "
      << std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
             .count()
      << "ms to complete CPU->GPU generation." << std::endl;
  m_built = true;

  mCRL2log(mcrl2::log::debug)
      << "Verts: " << m_total_vert << " Tris: " << m_total_tris << std::endl;
  mCRL2log(mcrl2::log::debug)
      << m_instances.size() << "/" << m_scenegraph.sceneData.nodes.size()
      << " draw calls" << std::endl;
}
/// @brief Initialize scene.
void Test::TScene::initializeScene()
{
  initializeOpenGLFunctions();
  glCheckError();
  int* viewport = new int[4];
  glGetIntegerv(GL_VIEWPORT, viewport);
  m_camera->setViewport(viewport[2] - viewport[0],
                        viewport[3] - viewport[1]);
  m_camera->reset();


  m_prog.addShaderFromSourceCode(QOpenGLShader::Vertex, VS);
  glCheckError();
  m_prog.addShaderFromSourceCode(QOpenGLShader::Fragment, FS);
  glCheckError();
  mCRL2log(mcrl2::log::debug) << "Program linked and bound." << std::endl;
  
  if (!m_prog.link())
  {
    mCRL2log(mcrl2::log::debug) << "Failed to link shader." << std::endl;
    mCRL2log(mcrl2::log::debug) << m_prog.log().toStdString() << std::endl;
  }

  if (!m_prog.bind())
  {
    mCRL2log(mcrl2::log::debug) << "Failed to bind shader." << std::endl;
    mCRL2log(mcrl2::log::debug) << m_prog.log().toStdString() << std::endl;
  }

  m_vao.create();
  m_vao.bind();

  m_vertex_loc = m_prog.attributeLocation("vertex");
  m_normal_loc = m_prog.attributeLocation("normal");
  m_color_loc = m_prog.attributeLocation("color");
  m_model_matrix_loc = m_prog.attributeLocation("model_mat");

  m_prog.enableAttributeArray(m_vertex_loc);
  m_prog.enableAttributeArray(m_normal_loc);
  m_prog.enableAttributeArray(m_color_loc);
  glVertexAttribDivisor(m_color_loc, 1);
  m_prog.enableAttributeArray(m_model_matrix_loc);
  for (int i = 0; i < 4; i++)
  {
    m_prog.enableAttributeArray(m_model_matrix_loc + i);
    glVertexAttribDivisor(m_model_matrix_loc + i, 1);
  }

  rebuildScene();
  // rebuildscene cleaned up
  mCRL2log(mcrl2::log::debug)
      << "[testscene] Device name: " << glGetString(GL_RENDERER)
      << std::endl;

  mCRL2log(mcrl2::log::debug)
      << "[testscene] " << m_prog.log().toStdString() << std::endl;

}

int color_offset = 0;
void Test::TScene::renderScene()
{
  //mCRL2log(mcrl2::log::debug) << "TScene::renderScene() " << m_built << std::endl;
  if (!m_built)
    rebuildScene();
  if (!m_prog.isLinked())
  {
    initializeScene();
  }

  m_prog.bind();
  glDisable(GL_BLEND);
  glCheckError();
  glDepthMask(GL_TRUE);
  glCheckError();
  glBlendFunc(GL_SRC_ALPHA, GL_ONE);
  glCheckError();
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LESS);
  glCheckError();
  glDisable(GL_CULL_FACE);
  glCheckError();
  //glCullFace(GL_BACK);
  glCheckError();

  QMatrix4x4 view = m_camera->getViewMatrix();
  QMatrix4x4 proj = m_camera->getProjectionMatrix();

#ifdef DEBUG_PRINT_VP
  mCRL2log(mcrl2::log::debug) << "View matrix: " << std::endl;
  for (int i = 0; i < 16; ++i)
  {
    mCRL2log(mcrl2::log::debug)
        << view.data()[i] << (((i + 1) % 4 == 0) ? "\n" : ", ");
  }
  mCRL2log(mcrl2::log::debug) << std::endl;

  mCRL2log(mcrl2::log::debug) << "Proj matrix: " << std::endl;
  for (int i = 0; i < 16; ++i)
  {
    mCRL2log(mcrl2::log::debug)
        << proj.data()[i] << (((i + 1) % 4 == 0) ? "\n" : ", ");
  }
  mCRL2log(mcrl2::log::debug) << std::endl;
#endif

  m_prog.setUniformValue("u_view", view);
  m_prog.setUniformValue("u_proj", proj);
  glCheckError();
  m_vao.bind();
  for (Instance* instance : m_instances)
  {
    instance->render(this, m_prog);
  }

  m_prog.release();
}