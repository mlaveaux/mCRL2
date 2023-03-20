// LTSView files
#include "scene.h"
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


const char* SSBO_VERTICES_NAME = "b_vertices";
const char* SSBO_NORMALS_NAME = "b_normals";
const char* SSBO_COLORS_NAME = "b_colors";
const char* SSBO_MATRICES_NAME = "b_matrices";

#ifndef NDEBUG
#define assertGL(msg) msgAssert(reportOpenGLError(), msg)
#else
#define assertGL(msg) noop
#endif

inline int cone_encode_resolution(int res, bool a, bool b)
{
  return res | (a ? (1 << 30) : 0) | (b ? (1 << 31) : 0);
}
int GlLTSView::ConeDB::getCone(GlUtil::Shapes::TruncatedCone* cone,
                               int resolution, SceneData& sceneData)
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
    GlLTSView::Mesh mesh =
        GlUtil::DefaultFactories::TruncatedConeFactory<Mesh>::createPrimitive(
            cone, resolution);
    int mesh_id = (int)sceneData.meshes.size();
    sceneData.meshes.push_back(mesh);
    cones[key] = mesh_id;
    return mesh_id;
  }
}

int GlLTSView::SphereDB::getSphere(GlUtil::Shapes::Sphere* sphere,
                                   int resolution, SceneData& sceneData)
{
  auto it = spheres.find({sphere->radius, resolution});
  if (it != spheres.end())
  {
    return it->second;
  }
  else
  {
    // create the mesh
    GlLTSView::Mesh mesh =
        GlUtil::DefaultFactories::SphereFactory<Mesh>::createPrimitive(
            sphere, resolution);
    int mesh_id = (int)sceneData.meshes.size();
    sceneData.meshes.push_back(mesh);
    spheres[{sphere->radius, resolution}] = mesh_id;
    return mesh_id;
  }
}

GlLTSView::Scene::Scene(QOpenGLWidget* glwidget, Cluster* root)
    : LTSScene(), m_clusterRoot(root)
{
}



/// TODO: Fix hardcoded resolutions
template <int sphereRes, int coneRes>
GlLTSView::SGNode* GlLTSView::SceneGraphFunctor<sphereRes, coneRes>::operator()(
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
        std::max((1+(coneRes/vistreenode->num_children))*vistreenode->num_children, 6*vistreenode->num_children), sceneData);
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

void GlLTSView::Scene::initializeScene()
{
  std::cout << "Initialize scene." << std::endl;
  // create shaderprogram
  if (!m_prog.addShaderFromSourceFile(
                   QOpenGLShader::Vertex,
                   "C:/Github/mCRL2/tools/release/ltsview/shaders/ltsview.vs"))
               std::cout << "Vertex shader can't be added." << std::endl;
  if (!m_prog.addShaderFromSourceFile(
                   QOpenGLShader::Fragment,
                   "C:/Github/mCRL2/tools/release/ltsview/shaders/ltsview.fs"))
               std::cout << "Fragment shader can't be added." << std::endl;

  
  if(!m_prog.link()){
        mCRL2log(mcrl2::log::error) << "Fatal error in linking program. Displaying log and exiting..." << std::endl;
        mCRL2log(mcrl2::log::error) << m_prog.log().toStdString() << std::endl;
        exit(-1);
    }
  if (!m_prog.bind()){
        mCRL2log(mcrl2::log::error) << "Fatal error in binding program. Displaying log and exiting..." << std::endl;
        mCRL2log(mcrl2::log::error) << m_prog.log().toStdString() << std::endl;
        exit(-1);
    }
  if (!built)
    rebuildScene();
  assertGL("Rebuild Scene");
  if (!built)
    return;
  QOpenGLExtraFunctions *f = QOpenGLContext::currentContext()->extraFunctions();


  GLint count;
  f->glGetProgramiv(m_prog.programId(), GL_ACTIVE_ATTRIBUTES, &count);
  GLint length, size;
  GLenum type;
  GLchar *name = (GLchar*) malloc(1024);
  mCRL2log(mcrl2::log::debug) << "Found: " << count << " attributes to be active." << std::endl;
  for (int i = 0; i < count; i++){
    f->glGetActiveAttrib(m_prog.programId(), (GLuint)i, 64, &length, &size, &type, name);
    mCRL2log(mcrl2::log::debug) << "Attr " << i << ": " << name << " length: " << length << " size: " << size << " type: " << type << std::endl;
  }


  // set all pointers to the right places
  m_uniform_view_loc = m_prog.uniformLocation("u_view");
  m_uniform_proj_loc = m_prog.uniformLocation("u_proj");
  m_uniform_alpha_loc = m_prog.uniformLocation("u_alpha");
  std::cout << "m_alpha_loc: " << m_uniform_alpha_loc << std::endl;
}

bool GlLTSView::Scene::reportOpenGLError()
{
  QOpenGLFunctions* f =
  QOpenGLContext::currentContext()->functions();
  GLenum error = f->glGetError();
  bool okay = error == GL_NO_ERROR;
  while (error != GL_NO_ERROR){
    if (error != GL_NO_ERROR)
    {
      std::string log = m_prog.log().toStdString();
      std::cout << "OpenGL error occured with code 0x" << std::hex << error
                << std::dec << (log.size() > 5 ? ":\n" + log : ".")
                << std::endl;
    }
    error = f->glGetError();
  }
  return okay;
}

static QVector4D color_to_vector(const QColor& color)
{
  return {
      static_cast<float>(color.redF()),
      static_cast<float>(color.greenF()),
      static_cast<float>(color.blueF()),
      static_cast<float>(color.alphaF())
         };
}

void GlLTSView::Scene::rebuildScene()
{
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
  if (built)
    free(m_scenegraph.root);
  if (built)
    free(m_vistreeRoot);
  if (built)
    m_scenegraph.sceneData = SceneData();
  int vertex_loc = m_prog.attributeLocation("vertex");
  int normal_loc = m_prog.attributeLocation("normal");
  int color_loc = m_prog.attributeLocation("color");
  int model_matrix_loc = m_prog.attributeLocation("model_mat");
  /// TODO: Change settings.h to have enum of modes instead of bool
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
                        GlLTSView::SceneGraphFunctor<sphereRes, coneRes>,
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
    new_instance.createModel(mesh, vertex_loc, normal_loc);
    
    n_verts += mesh.n_vertices;
    n_tris += mesh.n_triangles;
  }

  for (auto& node : m_scenegraph.sceneData.nodes)
  {
    Instance& instance = *m_instances[node.model_id];
    instance.color_cpu.push_back(color_to_vector(m_scenegraph.sceneData.colors[node.data_id]));
    instance.matrix_cpu.push_back(
        m_scenegraph.sceneData.matrices[node.data_id]);
    instance.instance_count++;
  }

  // enable attribute arrays
  m_prog.setAttributeBuffer(vertex_loc, GL_FLOAT, 0, 3);
  m_prog.enableAttributeArray(vertex_loc);

  m_prog.setAttributeBuffer(normal_loc, GL_FLOAT, 0, 3);
  m_prog.enableAttributeArray(normal_loc);

  m_prog.enableAttributeArray(color_loc);
  for (int i = 0; i < 4; i++)
  {
    m_prog.enableAttributeArray(model_matrix_loc + i);
  }

  for (auto& instance : m_instances)
  {
    instance->createInstances(color_loc, model_matrix_loc);
  }

  auto end = std::chrono::high_resolution_clock::now();

  mCRL2log(mcrl2::log::debug)
      << "Took "
      << std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
             .count()
      << "ms to complete CPU->GPU generation." << std::endl;
  built = true;
}

void GlLTSView::Scene::renderScene()
{
  QOpenGLExtraFunctions *f = QOpenGLContext::currentContext()->extraFunctions();
  if (!built)
    return;
  if (!m_scenegraph.root)
    return; /// TODO: Assert
  if (!buffers_exist)
    initializeScene();
    
  // QPainter _painter = QPainter(&m_glwidget);
  std::cout << "Render called" << std::endl;
  m_prog.bind();
  
  /// TODO: Check what needs to be updated.
  QMatrix4x4 view_matrix;
  view_matrix.setToIdentity();
  view_matrix.translate(-10, 0, 0);
  QMatrix4x4 proj_matrix;
  proj_matrix.perspective(45, 16.0 / 9.0, 0.1, 1000);
  m_prog.setUniformValue(m_uniform_view_loc, view_matrix);
  //glCheckError();
  assertGL("Set uniform; u_view");
  m_prog.setUniformValue(m_uniform_proj_loc, proj_matrix);
  //glCheckError();
  assertGL("Set uniform; u_proj");

  std::vector<const char*> names = {"b_vertices", "b_normals", "b_colors", "b_matrices", "u_view", "u_proj", "u_alpha", "in_vertexData"};
  std::vector<GLenum>      types = {GL_SHADER_STORAGE_BLOCK, GL_SHADER_STORAGE_BLOCK, GL_SHADER_STORAGE_BLOCK, GL_SHADER_STORAGE_BLOCK, GL_UNIFORM, GL_UNIFORM, GL_UNIFORM, GL_PROGRAM_INPUT};
  char* namebuff = new char[64];
  for (int i = 0 ; i < names.size(); i++){
    GLint resourceIndex = f->glGetProgramResourceIndex(m_prog.programId(), types[i], names[i]);
    GLint length;
    f->glGetProgramResourceName(m_prog.programId(), types[i], resourceIndex, 64, &length, namebuff);
    std::cout << "---- " << names[i] << " resource index: " << resourceIndex << " with resource name: " << namebuff << std::endl;
  }

  // set transparency
  m_prog.setUniformValue(m_uniform_alpha_loc,
                          Settings::instance().transparency.value() / 100.0f);

  assertGL("GlCullFace(GL_FRONT)");

  /// TODO: Multiple calls to get proper transparency.
  std::cout << "Before draw elements" << std::endl;

  for (auto& instance : m_instances)
  {
    instance->render();
  }

  m_prog.release();
  mCRL2log(mcrl2::log::debug) << "After draw elements" << std::endl;
  mCRL2log(mcrl2::log::debug) << "Painter done" << std::endl;
}