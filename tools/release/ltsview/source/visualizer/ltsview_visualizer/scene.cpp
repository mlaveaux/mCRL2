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
        std::max(coneRes, vistreenode->num_children), sceneData);
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
  if (!built)
    rebuildScene();
  if (!built)
    return;
  // create shaderprogram
  if (!program.addShaderFromSourceFile(
                   QOpenGLShader::Vertex,
                   "C:/Github/mCRL2/tools/release/ltsview/shaders/ltsview.vs"))
               std::cout << "Vertex shader can't be added." << std::endl;
  if (!program.addShaderFromSourceFile(
                   QOpenGLShader::Fragment,
                   "C:/Github/mCRL2/tools/release/ltsview/shaders/ltsview.fs"))
               std::cout << "Fragment shader can't be added." << std::endl;

  
  if (!program.link()) std::cout << "Program linking failed." << std::endl;
  if (!program.bind()) std::cout <<  "Program binding failed." << std::endl;;
  
  QOpenGLExtraFunctions *f = QOpenGLContext::currentContext()->extraFunctions();

  m_vao.create();
  QOpenGLVertexArrayObject::Binder vaoBinder(&m_vao);
  
  // generate the buffers with the currently active vao
  genBuffers();

  GLint count;
  f->glGetProgramiv(program.programId(), GL_ACTIVE_ATTRIBUTES, &count);
  GLint length, size;
  GLenum type;
  GLchar *name = (GLchar*) malloc(1024);
  std::cout << "Found: " << count << " attributes to be active." << std::endl;
  for (int i = 0; i < count; i++){
    f->glGetActiveAttrib(program.programId(), (GLuint)i, 64, &length, &size, &type, name);
    std::cout << "Attr " << i << ": " << name << " length: " << length << " size: " << size << " type: " << type << std::endl;
  }


  assertGL("Gen buffers function");

  // set all pointers to the right places
  m_view_loc = program.uniformLocation("u_view");
  m_proj_loc = program.uniformLocation("u_proj");
  m_alpha_loc = program.uniformLocation("u_alpha");
  std::cout << "m_alpha_loc: " << m_alpha_loc << std::endl;
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
      std::string log = program.log().toStdString();
      std::cout << "OpenGL error occured with code 0x" << std::hex << error
                << std::dec << (log.size() > 5 ? ":\n" + log : ".")
                << std::endl;
    }
    error = f->glGetError();
  }
  return okay;
}

void GlLTSView::Scene::createBufferObject(void* data, int num_bytes,
                                          QOpenGLBuffer& buff,
                                          const char* name, QOpenGLBuffer::UsagePattern usage_pattern)
{
  // make buffer
  buff.create();
  // bind it to vao
  buff.bind();
  buff.setUsagePattern(usage_pattern);
  
  // send data
  buff.allocate(data, num_bytes);

  QOpenGLExtraFunctions *f = QOpenGLContext::currentContext()->extraFunctions();

  GLuint loc = f->glGetProgramResourceIndex(program.programId(), GL_SHADER_STORAGE_BLOCK, name);

  // Bind buffer base
  f->glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, buff.bufferId());


  // Wait until writing is complete
  f->glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT); /// TODO: test if necessary


  /// TODO: Don't use << like this
  assertGL("Attempting to create buffer object: " << name);
}

void GlLTSView::Scene::genBuffers()
{
  QOpenGLExtraFunctions *f = QOpenGLContext::currentContext()->extraFunctions();
  std::cout << "Genbuffers called. " << std::endl;
  if (!built)
    return;
  assert(built); // we can't create buffers if there's no data
  if (buffers_exist)
  {
    buffers_exist = false;
  }
  int index = 0; // used for writing to correct places in arrays
  assert(!buffers_exist);

  // Create buffers
  // Fill VBO
  GLuint* vertexdata = new GLuint[2 * m_vertexData.size()];
  index = 0;
  for (auto& data : m_vertexData)
  {
    vertexdata[index++] = data.model_index;
    vertexdata[index++] = data.vertex_index;
  }
  m_vbo.create(); m_vbo.bind();
  m_vbo.setUsagePattern(QOpenGLBuffer::UsagePattern::StaticDraw);
  m_vbo.allocate(vertexdata, index*sizeof(GLuint));
  std::cout << "vbo enable vertexattribarray naar "
            << program.attributeLocation("in_vertexData") << " < QT | OpenGL > " << f->glGetAttribLocation(program.programId(), "in_vertexData") << std::endl;
  f->glEnableVertexAttribArray(0);
  f->glVertexAttribIPointer(0, 2, GL_UNSIGNED_INT, 0, 0);
  assertGL("VBO create");

  // Fill IBO
  GLuint* indices = new GLuint[3 * m_triangles.size()];
  index = 0;
  for (auto& tri : m_triangles)
  {
    indices[index++] = tri[0];
    indices[index++] = tri[1];
    indices[index++] = tri[2];
  }
  m_ibo = QOpenGLBuffer(QOpenGLBuffer::IndexBuffer);
  m_ibo.setUsagePattern(QOpenGLBuffer::UsagePattern::StaticDraw);
  m_ibo.allocate(indices, index*sizeof(GLuint));
  assertGL("IBO create");

  // Create vertex SSBO
  // temporarily copy data to gl array
  GLfloat* vertices =
      (GLfloat*)std::malloc(m_vertices.size() * 4 * sizeof(GLfloat));
  index = 0;
  for (auto& vert : m_vertices)
  {
    vertices[index++] = vert.x();
    vertices[index++] = vert.y();
    vertices[index++] = vert.z();
    vertices[index++] = 1;
  }
  createBufferObject(vertices, index * sizeof(GLfloat), m_vertex_ssbo,
                     SSBO_VERTICES_NAME, QOpenGLBuffer::UsagePattern::StaticDraw);
  glCheckError();

  // Create normal SSBO
  // temporarily copy data to gl array
  GLfloat* normals =
      (GLfloat*)std::malloc(m_normals.size() * 4 * sizeof(GLfloat));
  index = 0;
  for (auto& norm : m_normals)
  {
    normals[index++] = norm.x();
    normals[index++] = norm.y();
    normals[index++] = norm.z();
    normals[index++] = 0; // important for transforming normals
  }
  createBufferObject(normals, index * sizeof(GLfloat), m_normal_ssbo,
                     SSBO_NORMALS_NAME, QOpenGLBuffer::UsagePattern::StaticDraw);
  glCheckError();

  // Create color SSBO
  // temporarily copy data to gl array
  GLfloat* colors = (GLfloat*)std::malloc(m_scenegraph.sceneData.colors.size() *
                                          4 * sizeof(GLfloat));
  index = 0;
  for (auto& color : m_scenegraph.sceneData.colors)
  {
    colors[index++] = color.redF();
    colors[index++] = color.greenF();
    colors[index++] = color.blueF();
    colors[index++] = 1;
  }
  // NB: GL_DYNAMIC_DRAW because colors can be updated at will
  createBufferObject(colors, index * sizeof(GLfloat), m_color_ssbo,
                     SSBO_COLORS_NAME, QOpenGLBuffer::UsagePattern::StaticDraw);
  glCheckError();

  // Create matrix SSBO
  GLfloat* matrices = (GLfloat*)std::malloc(
      m_scenegraph.sceneData.matrices.size() * 16 * sizeof(GLfloat));
  GLfloat* loc = matrices;
  for (auto& mat : m_scenegraph.sceneData.matrices)
  {
    mat.copyDataTo(loc);
    loc += 16;
  }
  // NB: Like with colors, can be updated at will so dynamic draw hint
  createBufferObject(
      matrices, m_scenegraph.sceneData.matrices.size() * 16 * sizeof(GLfloat),
      m_matrix_ssbo, SSBO_MATRICES_NAME, QOpenGLBuffer::UsagePattern::StaticDraw);
  glCheckError();

  assertGL("Create all buffers.");
  buffers_exist = true;

  // Delete all temporary data arrays
  delete vertexdata;
  delete indices;
  delete vertices;
  delete matrices;
  delete normals;
  delete colors;
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
  if (built)
    free(m_scenegraph.root);
  if (built)
    free(m_vistreeRoot);
  if (built)
    m_scenegraph.sceneData = SceneData();
  /// TODO: Change settings.h to have enum of modes instead of bool
  m_vistreeRoot = VisTreeGenerator::generate(
      Settings::instance().clusterVisStyleTubes.value()
          ? VisTreeGenerator::Mode::TUBES
          : VisTreeGenerator::Mode::CONES,
      m_clusterRoot);
  const int sphereRes = 3;
  const int coneRes = 20;
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
  for (const NodeData& node : m_scenegraph.sceneData.nodes)
  {
    n_verts += m_scenegraph.sceneData.meshes[node.model_id].n_vertices;
    n_tris += m_scenegraph.sceneData.meshes[node.model_id].n_triangles;
  }

  m_vertices = std::vector<QVector3D>();
  m_vertices.reserve(n_verts);
  auto back_inserter_vertices = std::back_inserter(m_vertices);

  m_normals = std::vector<QVector3D>();
  m_normals.reserve(n_verts);
  auto back_inserter_normals = std::back_inserter(m_normals);

  m_triangles = std::vector<GlUtil::Triangle>();
  m_triangles.reserve(n_tris);
  auto back_inserter_triangles = std::back_inserter(m_triangles);

  std::cout << "Nodes: " << m_scenegraph.sceneData.nodes.size() << std::endl;
  std::cout << "Meshes: " << m_scenegraph.sceneData.meshes.size() << std::endl;

  int compr_verts = 0, compr_tris = 0;
  for (auto& mesh : m_scenegraph.sceneData.meshes)
  {
    compr_verts += mesh.n_vertices;
    compr_tris += mesh.n_triangles;
  }
  std::cout << "Stored verts: " << compr_verts
            << ". Stored tris: " << compr_tris << std::endl;
  std::cout << "Total verts: " << n_verts << ". Total tris: " << n_tris
            << std::endl;

  int vert_offset = 0;
  int tri_offset = 0;
  std::vector<QColor> colors = std::vector<QColor>(n_verts);
  for (auto& mesh : m_scenegraph.sceneData.meshes)
  {
    assert(tri_offset < n_tris);
    assert(vert_offset < n_verts);

    std::move(mesh.vertices.begin(), mesh.vertices.end(),
              back_inserter_vertices);
    std::move(mesh.vertex_normals.begin(), mesh.vertex_normals.end(),
              back_inserter_normals);
    std::move(mesh.triangles.begin(), mesh.triangles.end(),
              back_inserter_triangles);

    for (int i = tri_offset; i < tri_offset + mesh.n_triangles; ++i)
    {
      m_triangles[i][0] += vert_offset;
      m_triangles[i][1] += vert_offset;
      m_triangles[i][2] += vert_offset;

      assert(m_triangles[i][0] <= n_verts);
      assert(m_triangles[i][1] <= n_verts);
      assert(m_triangles[i][2] <= n_verts);
    }
    tri_offset += mesh.n_triangles;
    vert_offset += mesh.n_vertices;
  }
  auto end = std::chrono::high_resolution_clock::now();

  std::cout << "Took "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                     start)
                   .count()
            << "ms to complete." << std::endl;
  long uncompr_bytes = n_verts * 96;
  std::cout << "Used resolutions: " << std::endl;
  std::cout << "\t- Sphere: " << sphereRes << std::endl;
  std::cout << "\t- Cone: " << coneRes << " (adaptive)" << std::endl;
  std::cout << uncompr_bytes
            << " bytes in total for vertices only with no compression with "
               "pos/color/transform: "
            << std::endl;
  std::cout << "\t- " << uncompr_bytes / 1024.0 << "KB" << std::endl;
  std::cout << "\t- " << uncompr_bytes / (1024.0 * 1024.0) << "MB" << std::endl;
  std::cout << "\t- " << uncompr_bytes / (1024.0 * 1024.0 * 1024) << "GB"
            << std::endl;
  long compr_bytes =
      compr_verts * 16 + m_scenegraph.sceneData.meshes.size() * (64 + 16);
  std::cout << compr_bytes
            << " bytes in total with compressed storage: " << std::endl;
  std::cout << "\t- " << compr_bytes / 1024.0 << "KB" << std::endl;
  std::cout << "\t- " << compr_bytes / (1024.0 * 1024.0) << "MB" << std::endl;
  std::cout << "\t- " << compr_bytes / (1024.0 * 1024.0 * 1024) << "GB"
            << std::endl;

  auto sphereMesh = GlUtil::DefaultFactories::
      SphereFactory<GlUtil::Meshes::MeshTypes::TriangleMesh>::createPrimitive(
          new GlUtil::Shapes::Sphere(), sphereRes);
  sphereMesh.toObj("sphere-res=" + std::to_string(sphereRes) + ".obj");
  auto cone = new GlUtil::Shapes::TruncatedCone();
  cone->radius_top = 1;
  cone->radius_bot = 2;
  cone->height = 1;
  cone->fill_bot = true;
  cone->fill_top = true;
  auto coneMesh = GlUtil::DefaultFactories::TruncatedConeFactory<
      GlUtil::Meshes::MeshTypes::TriangleMesh>::createPrimitive(cone, coneRes);
  coneMesh.toObj("cone-res=" + std::to_string(coneRes) + ".obj");
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
  program.bind();
  
  /// TODO: Check what needs to be updated.
  QMatrix4x4 view_matrix;
  view_matrix.setToIdentity();
  view_matrix.translate(-10, 0, 0);
  QMatrix4x4 proj_matrix;
  proj_matrix.perspective(45, 16.0 / 9.0, 0.1, 1000);
  program.setUniformValue(m_view_loc, view_matrix);
  glCheckError();
  assertGL("Set uniform; u_view");
  program.setUniformValue(m_proj_loc, proj_matrix);
  glCheckError();
  assertGL("Set uniform; u_proj");

  std::vector<const char*> names = {"b_vertices", "b_normals", "b_colors", "b_matrices", "u_view", "u_proj", "u_alpha", "in_vertexData"};
  std::vector<GLenum>      types = {GL_SHADER_STORAGE_BLOCK, GL_SHADER_STORAGE_BLOCK, GL_SHADER_STORAGE_BLOCK, GL_SHADER_STORAGE_BLOCK, GL_UNIFORM, GL_UNIFORM, GL_UNIFORM, GL_PROGRAM_INPUT};
  char* namebuff = new char[64];
  for (int i = 0 ; i < names.size(); i++){
    GLint resourceIndex = f->glGetProgramResourceIndex(program.programId(), types[i], names[i]);
    GLint length;
    f->glGetProgramResourceName(program.programId(), types[i], resourceIndex, 64, &length, namebuff);
    std::cout << "---- " << names[i] << " resource index: " << resourceIndex << " with resource name: " << namebuff << std::endl;
  }

  // set transparency
  program.setUniformValue(m_alpha_loc,
                          Settings::instance().transparency.value() / 100.0f);

  f->glClearColor(Settings::instance().backgroundColor.value().red() / 255.0,
               Settings::instance().backgroundColor.value().green() / 255.0,
               Settings::instance().backgroundColor.value().blue() / 255.0,
               1.0f);
  // Clear existing buffers
  f->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  f->glEnable(GL_DEPTH_TEST);
  f->glDepthFunc(GL_LESS); // closer to camera wins

  f->glCullFace(GL_BACK);
  assertGL("GlCullFace(GL_FRONT)");

  /// TODO: Multiple calls to get proper transparency.
  std::cout << "Before draw elements" << std::endl;

  QOpenGLVertexArrayObject::Binder vaoBinder(&m_vao); // bind vao
  std::cout << "bind vao: " << m_vao.isCreated() << std::endl;
  glCheckError();
  GLint count;
  f->glGetProgramiv(program.programId(), GL_ACTIVE_ATTRIBUTES, &count);
  GLint length, size;
  GLenum type;
  GLchar *name = (GLchar*) malloc(64);
  std::cout << "Found: " << count << " attributes to be active." << std::endl;
  for (int i = 0; i < count; i++){
    f->glGetActiveAttrib(program.programId(), (GLuint)i, 64, &length, &size, &type, name);
    
    std::cout << "Attr " << i << ": " << name << " length: " << length << " size: " << size << " type: " << (type == GL_INT_VEC2 ? "GL_INT_VEC2" : "?????") << std::endl;
  }
  std::cout << "m_vbo: " << m_vbo.bufferId() << std::endl;
  std::cout << "m_ibo: " << m_ibo.bufferId() << std::endl;
  f->glDrawElements(GL_TRIANGLES, 3, GL_UNSIGNED_INT, (void*)0);

  glCheckError();
  std::cout << "After draw elements" << std::endl;
  glCheckError();
  std::cout << "Painter done" << std::endl;
  exit(0);
}