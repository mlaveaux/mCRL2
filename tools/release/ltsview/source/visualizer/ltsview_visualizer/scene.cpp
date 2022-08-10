// LTSView files
#include "scene.h"
#include "cluster.h"
#include "vistreegenerator.h"
#include "settings.h"
// Generic opengl files
#include "glutil.h"
#include "glscenegraph.h"
#include "glprimitivefactories.h"



GlLTSView::Scene::Scene(QOpenGLWidget& glwidget, SceneGraph<GlLTSView::NodeData, GlLTSView::SceneData>& scenegraph, Camera& camera, Cluster *root) : LTSScene(glwidget, scenegraph, camera), m_clusterRoot(root){}


/// TODO: Fix hardcoded resolutions
/// TODO: Use conedb/spheredb
GlLTSView::SGNode* GlLTSView::SceneGraphFunctor::operator()(SGNode* parent, VisTree::VisTreeNode* vistreenode){
    SGNode* node = new SGNode();
    Mesh mesh;
    QMatrix4x4 matrix;
    switch (vistreenode->data.shape->getShapeType()){
        case GlUtil::ShapeType::SPHERE:
      std::cout << "building sphere" << std::endl;
            mesh = GlUtil::DefaultFactories::SphereFactory<Mesh>::createPrimitive(static_cast<GlUtil::Shapes::Sphere*>(vistreenode->data.shape), 2);
            matrix = vistreenode->data.matrix;
            node->data.model_id = sceneData.meshes.size();
            node->data.model_matrix_id = sceneData.matrices.size();
            sceneData.meshes.emplace_back(mesh);
            sceneData.matrices.emplace_back(matrix);
            return node;
        case GlUtil::ShapeType::TRUNCATED_CONE:
          std::cout << "building cone" << std::endl;
            mesh = GlUtil::DefaultFactories::TruncatedConeFactory<Mesh>::createPrimitive(static_cast<GlUtil::Shapes::TruncatedCone*>(vistreenode->data.shape), 20);
            matrix = vistreenode->data.matrix;
            node->data.model_id = sceneData.meshes.size();
            node->data.model_matrix_id = sceneData.matrices.size();
            sceneData.meshes.emplace_back(mesh);
            sceneData.matrices.emplace_back(matrix);
            return node;
        default:
            /// TODO: Log unknown shape
            node->data.model_id = -1;
            node->data.model_matrix_id = -1;
            return node;
    }
    return node; 
}

void GlLTSView::Scene::rebuild(){
    // clear the scene data
  std::cout << "Rebuild called." << std::endl;
    m_scenegraph.sceneData.meshes.clear(); m_scenegraph.sceneData.matrices.clear();
    if (built) free(m_scenegraph.root);
    if (built) free(m_vistreeRoot);

    /// TODO: Change settings.h to have enum of modes instead of bool
    m_vistreeRoot = VisTreeGenerator::generate(Settings::instance().clusterVisStyleTubes.value() ? VisTreeGenerator::Mode::TUBES : VisTreeGenerator::Mode::CONES, m_clusterRoot);

    std::cout << "Creating functor." << std::endl;
    SceneGraphFunctor sgf(m_scenegraph.sceneData);

    m_scenegraph.root = GlUtil::fold_tree<SGNode, VisTree::VisTreeNode, GlLTSView::SceneGraphFunctor, VisTree::VisTreeNode::childIterator>(
        nullptr, m_vistreeRoot,
        SGNode::setParent, SGNode::addChild,
        sgf, VisTree::VisTreeNode::getChildBegin, VisTree::VisTreeNode::getChildEnd);

    int n_verts = 0;
    int n_tris  = 0; 
    for (const Mesh& mesh : m_scenegraph.sceneData.meshes){
        n_verts += mesh.n_vertices;
        n_tris  += mesh.n_triangles;
    } 

    std::cout << "Vertices: " << n_verts << " tris: " << n_tris << std::endl;


    program.addShaderFromSourceFile(QOpenGLShader::Vertex, "C:/Github/mCRL2/tools/release/ltsview/source/generic_visualizer/shaders/simple.vs");
    program.addShaderFromSourceFile(QOpenGLShader::Fragment, "C:/Github/mCRL2/tools/release/ltsview/source/generic_visualizer/shaders/simple.fs");

    program.link();
    program.bind();


    // Generate vertex arrays and bind
    glGenVertexArrays(1, &m_VAO);
    glBindVertexArray(m_VAO);

    int total_vertex_data = std::accumulate(m_scenegraph.sceneData.meshes.begin(), m_scenegraph.sceneData.meshes.end(), 0, [](int accum, auto& mesh) { return accum + mesh.n_vertices; });
    std::vector<QVector3D> vertices(total_vertex_data);
    auto back_inserter_vertices = std::back_inserter(vertices);
    std::vector<QVector3D> normals(total_vertex_data);
    auto back_inserter_normals = std::back_inserter(vertices);
    for(auto& mesh : m_scenegraph.sceneData.meshes){
        std::move(mesh.vertices.begin(), mesh.vertices.end(), back_inserter_vertices);
        std::move(mesh.vertex_normals.begin(), mesh.vertex_normals.end(), back_inserter_normals);
    } 
    
    total_vertex_data *= 3 * sizeof(float);
    glGenBuffers(1, &m_PositionBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, m_PositionBuffer);
    glBufferData(GL_ARRAY_BUFFER, total_vertex_data, &vertices[0], GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(0);

    glGenBuffers(1, &m_NormalBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, m_NormalBuffer);
    glBufferData(GL_ARRAY_BUFFER, total_vertex_data, &normals[0], GL_STATIC_DRAW);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_TRUE, 0, 0);
    glEnableVertexAttribArray(1);

    built = true;
}

void GlLTSView::Scene::render(QPainter& painter){
//    m_renderer.render(painter);
}