#include "testscene.h"
#include "cluster.h"
#include "vistreegenerator.h"
#include "settings.h"
// Generic opengl files
#include "glutil.h"
#include "glscenegraph.h"
#include "glprimitivefactories.h"

// std includes
#include <fstream>
#include <sstream>
#include <chrono>

#include <QOpenGLFunctions>
#include <QOpenGLExtraFunctions>
#include <QOpenGLFunctions_4_3_Core>

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
        std::max((1+(coneRes/vistreenode->num_children))*vistreenode->num_children, 6*vistreenode->num_children)+1, sceneData);
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
    "#version 430 compatibility\n"

    "layout(std430) buffer b_verts { vec4 vertices[]; };\n"

    "layout(std430) buffer b_normals { vec4 normals[]; };\n"

    "layout(std430) buffer b_colors { vec4 colors[]; };\n"

    "layout(std430) buffer b_matrices { mat4 matrices[]; };\n"
    // vert.x = transform/color
    // vert.y = vertex index
    "layout(location = 0) in ivec2 vert;\n"
    "out vec4 vColor;\n"

    "uniform mat4 u_view;\n"

    "uniform mat4 u_proj;\n"

    "out vec4 vNormal;\n"

    "void main(void)\n"
    "{\n"
    "   gl_Position = u_proj*u_view*matrices[vert.x]*vertices[vert.y];\n"
    "   vColor      =  colors[vert.x];\n"
    "   vNormal     = normals[vert.x];\n"
    "}";

const char* FS = "#version 430 compatibility\n"
                 "in vec4 vNormal;"
                 "in vec4 vColor;\n"
                 "out vec4 fColor;\n"
                 "void main(void)\n"
                 "{\n"
                 "   fColor = vColor;\n"
                 "}";

Test::TScene::TScene(Cluster* root) : TestScene()
{
  m_clusterRoot = root;
}

void Test::TScene::createBufferObject(QOpenGLBuffer& buff,
                                      const char* buff_name, const void* data,
                                      std::size_t size,
                                      QOpenGLExtraFunctions* f,
                                      QOpenGLFunctions_4_3_Core* f430)
{
  buff.create();
  buff.bind();
  buff.allocate(data, size);
  mCRL2log(mcrl2::log::debug)
      << "Buffer \"" << buff_name << "\""
      << (buff.isCreated() ? " created succesfully" : " failed to create")
      << "; buffer id : " << buff.bufferId() << std::endl;

  GLuint buffer_index = f->glGetProgramResourceIndex(
      m_prog.programId(), GL_SHADER_STORAGE_BLOCK, buff_name);
  f430->glShaderStorageBlockBinding(m_prog.programId(), buffer_index,
                                    m_ssbo_count);
  f->glBindBufferBase(GL_SHADER_STORAGE_BUFFER, m_ssbo_count++,
                      buff.bufferId());
  mCRL2log(mcrl2::log::debug)
      << "Bound \"" << buff_name << "\" to buffer_index: " << buffer_index
      << std::endl;
  buff.release();
}

void Test::TScene::reallocateStorage()
{
  if (m_cpu_allocated)
  {
    free(m_vertices_cpu);
    free(m_normals_cpu);
    free(m_triangles_cpu);
    free(m_vertexdata_cpu);
    free(m_matrices_cpu);
    free(m_colors_cpu);

    m_prog.removeAllShaders();

    m_cpu_allocated = false;
  }

  if (m_total_vert >= 0 && m_total_tris >= 0)
  {

    m_vertices_cpu = (float*)std::malloc(m_total_vert * sizeof(float) * 4);

    m_normals_cpu = (float*)std::malloc(m_total_vert * sizeof(float) * 4);

    m_triangles_cpu = (int*)std::malloc(m_total_tris * sizeof(int) * 3);

    m_vertexdata_cpu = (int*)std::malloc(m_total_vert * sizeof(int) * 2);

    m_matrices_cpu = (float*)std::malloc(m_scenegraph.sceneData.nodes.size() *
                                         sizeof(float) * 16);
    m_colors_cpu = (float*)std::malloc(m_scenegraph.sceneData.nodes.size() *
                                       sizeof(float) * 4);

    m_cpu_allocated = true;
  }
}

void Test::TScene::rebuildScene()
{
  if (m_vao.isCreated() && m_vao.objectId() != -1)
    m_vao.destroy();
  if (m_vbo.isCreated() && m_vbo.bufferId() != -1)
    m_vbo.destroy();
  if (m_ibo.isCreated() && m_ibo.bufferId() != -1)
    m_ibo.destroy();
  if (m_ssbo1.isCreated() && m_ssbo1.bufferId() != -1)
    m_ssbo1.destroy();
  if (m_ssbo2.isCreated() && m_ssbo2.bufferId() != -1)
    m_ssbo2.destroy();
  if (m_ssbo3.isCreated() && m_ssbo3.bufferId() != -1)
    m_ssbo3.destroy();
  if (m_ssbo4.isCreated() && m_ssbo4.bufferId() != -1)
    m_ssbo4.destroy();
  m_ssbo_count = 0;
  if (m_camera)
    m_camera->reset();
  if (!m_clusterRoot)
    return;
  assert(m_clusterRoot);

  m_scenegraph.sceneData.clear();


  //if (m_scenegraph.root)
    //delete m_scenegraph.root;

  /// TODO: Verify delete[] vs delete with someone who has more knowledge - Ruben
  //if (m_vistreeRoot)
    //delete m_vistreeRoot;

  /// TODO: Change settings.h to have enum of modes instead of bool
  m_vistreeRoot = VisTreeGenerator::generate(
      Settings::instance().clusterVisStyleTubes.value()
          ? VisTreeGenerator::Mode::TUBES
          : VisTreeGenerator::Mode::CONES,
      m_clusterRoot);

  /// TODO: Expose and/or make user configurable
  const int sphereRes = 3;
  const int coneRes = 1000;
  SceneGraphFunctor<sphereRes, coneRes> sgf(m_scenegraph.sceneData);

  m_scenegraph.root =
      GlUtil::fold_tree<SGNode, VisTree::VisTreeNode,
                        Test::SceneGraphFunctor<sphereRes, coneRes>,
                        VisTree::VisTreeNode::childIterator>(
          nullptr, m_vistreeRoot, SGNode::setParent, SGNode::addChild, sgf,
          VisTree::VisTreeNode::getChildBegin,
          VisTree::VisTreeNode::getChildEnd);

  // Scenegraph is now built as a vistree
  m_total_vert = 0;
  m_total_tris = 0;
  for (const NodeData& node : m_scenegraph.sceneData.nodes)
  {
    m_total_vert += m_scenegraph.sceneData.meshes[node.model_id].n_vertices;
    m_total_tris += m_scenegraph.sceneData.meshes[node.model_id].n_triangles;
  }

  if (!(m_total_vert > 0 && m_total_tris > 0))
  {
    // cant fill storage with nothing
    return;
  }

  reallocateStorage();

  std::vector<int> vertex_offsets({0});
  for (const auto& mesh : m_scenegraph.sceneData.meshes)
  {
    vertex_offsets.push_back(vertex_offsets.back() + mesh.n_vertices);
  }

  int vertexdata_index = 0;
  int tri_index = 0;
  for (const NodeData& node : m_scenegraph.sceneData.nodes)
  {
    Mesh& mesh = m_scenegraph.sceneData.meshes[node.model_id];
    int vert_offset = vertex_offsets[node.model_id];
    int tri_offset = vertexdata_index / 2;
    for (int i = 0; i < mesh.n_vertices; ++i)
    {
      m_vertexdata_cpu[vertexdata_index++] = node.data_id;
      m_vertexdata_cpu[vertexdata_index++] = vert_offset + i;
    }

    for (const auto& tri : mesh.triangles)
    {
      m_triangles_cpu[tri_index++] = tri[0] + tri_offset;
      m_triangles_cpu[tri_index++] = tri[1] + tri_offset;
      m_triangles_cpu[tri_index++] = tri[2] + tri_offset;
    }
  }

  int vert_index = 0;
  int norm_index = 0;

  for (auto& mesh : m_scenegraph.sceneData.meshes)
  {
    for (const QVector3D& vec : mesh.vertices)
    {
      m_vertices_cpu[vert_index++] = vec.x();
      m_vertices_cpu[vert_index++] = vec.y();
      m_vertices_cpu[vert_index++] = vec.z();
      m_vertices_cpu[vert_index++] = 1; // important to have w = 1 for positions
#ifdef DEBUG_PRINT_VERTICES
      mCRL2log(mcrl2::log::debug)
          << "(" << vec.x() << ", " << vec.y() << ", " << vec.z() << ", "
          << "1)" << std::endl;
#endif
    }
    for (const QVector3D& vec : mesh.vertex_normals)
    {
      m_normals_cpu[norm_index++] = vec.x();
      m_normals_cpu[norm_index++] = vec.y();
      m_normals_cpu[norm_index++] = vec.z();
      m_normals_cpu[norm_index++] = 0; // important to have w = 0 for normals
    }
  }

  int mat_index = 0;
  for (auto& mat : m_scenegraph.sceneData.matrices)
  {
    std::memcpy(m_matrices_cpu + mat_index, mat.data(), 16 * sizeof(float));
#ifdef DEBUG_PRINT_MATRICES
    mCRL2log(mcrl2::log::debug) << "Matrix: " << std::endl;
    for (int i = mat_index; i < mat_index + 16; ++i)
    {
      mCRL2log(mcrl2::log::debug)
          << m_matrices_cpu[i] << (((i + 1) % 4) != 0 ? ", " : "\n");
    }
#endif
    mat_index += 16;
  }

  int col_index = 0;
  for (auto& col : m_scenegraph.sceneData.colors)
  {
    m_colors_cpu[col_index++] = col.redF();
    m_colors_cpu[col_index++] = col.blueF();
    m_colors_cpu[col_index++] = col.greenF();
    m_colors_cpu[col_index++] = col.alphaF();
  }
  m_built = true;
}
/// @brief Initialize scene. Assumption that initializeOpenGLFunctions() has
/// been called in advance.
void Test::TScene::initializeScene()
{
  rebuildScene();
  if (!(m_total_vert > 0 && m_total_tris > 0))
    return;
  else
  {
    mCRL2log(mcrl2::log::debug)
        << "Verts: " << m_total_vert << " Tris: " << m_total_tris << std::endl;
  }
  // rebuildscene cleaned up

  m_prog.addShaderFromSourceCode(QOpenGLShader::Vertex, VS);
  m_prog.addShaderFromSourceCode(QOpenGLShader::Fragment, FS);
  m_prog.link();

  m_prog.bind();

  mCRL2log(mcrl2::log::debug) << "Program linked and bound." << std::endl;

  QOpenGLExtraFunctions* f = QOpenGLContext::currentContext()->extraFunctions();
  QOpenGLFunctions_4_3_Core* f430 =
      QOpenGLContext::currentContext()
          ->versionFunctions<QOpenGLFunctions_4_3_Core>();
  m_vao.create();
  QOpenGLVertexArrayObject::Binder vaoBinder(&m_vao);

  // VBO
  m_vbo.create();
  m_vbo.bind();
  m_vbo.allocate(m_vertexdata_cpu, 2 * sizeof(int) * m_total_vert);
  f->glEnableVertexAttribArray(0);
  f->glVertexAttribIPointer(0, 2, GL_INT, 0, 0);
  m_vbo.release();

  // IBO
  m_ibo = QOpenGLBuffer(QOpenGLBuffer::IndexBuffer);
  m_ibo.create();
  m_ibo.bind();
  m_ibo.allocate(m_triangles_cpu, 3 * sizeof(int) * m_total_tris);
  m_ibo.release();

  // SSBO
  createBufferObject(m_ssbo1, "b_verts", m_vertices_cpu,
                     m_total_vert * 4 * sizeof(float), f, f430);

  mCRL2log(mcrl2::log::debug) << "SSBO1: " << m_ssbo1.isCreated() << ", "
                              << m_ssbo1.bufferId() << std::endl;

  createBufferObject(m_ssbo2, "b_normals", m_normals_cpu,
                     m_total_vert * 4 * sizeof(float), f, f430);

  mCRL2log(mcrl2::log::debug) << "SSBO2: " << m_ssbo2.isCreated() << ", "
                              << m_ssbo2.bufferId() << std::endl;

  createBufferObject(m_ssbo3, "b_colors", m_colors_cpu,
                     4 * m_scenegraph.sceneData.colors.size() * sizeof(float),
                     f, f430);

  mCRL2log(mcrl2::log::debug) << "SSBO3: " << m_ssbo3.isCreated() << ", "
                              << m_ssbo3.bufferId() << std::endl;

  createBufferObject(
      m_ssbo4, "b_matrices", m_matrices_cpu,
      16 * m_scenegraph.sceneData.matrices.size() * sizeof(float), f,
      f430); // column major

  mCRL2log(mcrl2::log::debug) << "SSBO4: " << m_ssbo4.isCreated() << ", "
                              << m_ssbo4.bufferId() << std::endl;

  glCheckError();

  mCRL2log(mcrl2::log::debug)
      << "[testscene] Device name: " << f->glGetString(GL_RENDERER)
      << std::endl;

  mCRL2log(mcrl2::log::debug)
      << "[testscene] " << m_prog.log().toStdString() << std::endl;
}

int color_offset = 0;
void Test::TScene::renderScene()
{
  if (!m_built)
    return;
  if (!m_prog.isLinked())
  {
    initializeScene();
  }
  QOpenGLExtraFunctions* f = QOpenGLContext::currentContext()->extraFunctions();
  f->glEnable(GL_BLEND);
  glCheckError();
  f->glDepthMask(GL_TRUE);
  glCheckError();
  f->glBlendFunc(GL_SRC_ALPHA, GL_ONE);
  glCheckError();
  f->glDisable(GL_DEPTH_TEST);
  glCheckError();
  f->glEnable(GL_CULL_FACE);
  glCheckError();
  f->glCullFace(GL_BACK);
  glCheckError();

  m_prog.bind();

  QMatrix4x4 view = m_camera->getViewMatrix();
  QMatrix4x4 proj = m_camera->getProjectionMatrix();

#ifdef DEBUG_PRINT_PV
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

  QOpenGLVertexArrayObject::Binder vaoBinder(&m_vao);

  m_ibo.bind();
  f->glDrawElements(GL_TRIANGLES, m_total_tris * 3, GL_UNSIGNED_INT, (void*)0);
  glCheckError();
  vaoBinder.release();
  m_prog.release();
}