// LTSView files
#include "scene.h"
#include "cluster.h"
#include "vistreegenerator.h"
#include "settings.h"
// Generic opengl files
#include "glutil.h"
#include "glscenegraph.h"
#include "glprimitivefactories.h"

#include <fstream>



GlLTSView::Scene::Scene(QOpenGLWidget& glwidget, SceneGraph<GlLTSView::NodeData, GlLTSView::SceneData>& scenegraph, Camera& camera, Cluster *root) : LTSScene(glwidget, scenegraph, camera), m_clusterRoot(root){
    // initializeOpenGLFunctions();
}


/// TODO: Fix hardcoded resolutions
/// TODO: Use conedb/spheredb
GlLTSView::SGNode* GlLTSView::SceneGraphFunctor::operator()(SGNode* parent, VisTree::VisTreeNode* vistreenode){
    SGNode* node = new SGNode();
    Mesh mesh;
    QMatrix4x4 matrix;
    switch (vistreenode->data.shape->getShapeType()){
        case GlUtil::ShapeType::SPHERE:
            mesh = GlUtil::DefaultFactories::SphereFactory<Mesh>::createPrimitive(static_cast<GlUtil::Shapes::Sphere*>(vistreenode->data.shape), 3);
            matrix = vistreenode->data.matrix;
            node->data.model_id = sceneData.meshes.size();
            node->data.model_matrix_id = sceneData.matrices.size();
            sceneData.meshes.emplace_back(mesh);
            sceneData.matrices.emplace_back(matrix);
            return node;
        case GlUtil::ShapeType::TRUNCATED_CONE:
            mesh = GlUtil::DefaultFactories::TruncatedConeFactory<Mesh>::createPrimitive(static_cast<GlUtil::Shapes::TruncatedCone*>(vistreenode->data.shape), 50);
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

    std::vector<QVector3D> vertices = std::vector<QVector3D>();
    vertices.reserve(n_verts);
    auto back_inserter_vertices = std::back_inserter(vertices);
    std::vector<QVector3D> normals = std::vector<QVector3D>();
    normals.reserve(n_verts);
    auto back_inserter_normals = std::back_inserter(normals);
    std::vector<GlUtil::Triangle> triangles = std::vector<GlUtil::Triangle>();
    triangles.reserve(n_tris);

    auto back_inserter_triangles = std::back_inserter(triangles);
    int vert_offset = 0;
    int tri_offset = 0;
    for(int j = 0; j < m_scenegraph.sceneData.meshes.size(); ++j){
        auto& mesh = m_scenegraph.sceneData.meshes[j];
        auto& matrix = m_scenegraph.sceneData.matrices[j];
        std::move(mesh.vertices.begin(), mesh.vertices.end(), back_inserter_vertices);
        std::move(mesh.vertex_normals.begin(), mesh.vertex_normals.end(), back_inserter_normals);
        for (int i = vert_offset; i < vert_offset + mesh.n_vertices; i++){
            vertices[i] = matrix * vertices[i];
        }
        
        std::move(mesh.triangles.begin(), mesh.triangles.end(), back_inserter_triangles);
        for (int i = tri_offset; i < tri_offset + mesh.n_triangles; ++i){
            triangles[i][0] += vert_offset;
            triangles[i][1] += vert_offset;
            triangles[i][2] += vert_offset;

            if (triangles[i][0] >= n_verts) std::cout << "Problem with tris. " << std::endl;
            if (triangles[i][1] >= n_verts) std::cout << "Problem with tris. " << std::endl;
            if (triangles[i][2] >= n_verts) std::cout << "Problem with tris. " << std::endl;
        }
        tri_offset += mesh.n_triangles;
        vert_offset += mesh.n_vertices;
    }

    std::ofstream file("output.obj");
    if (!file)
    {
        return ;
    }


    std::cout << "Vertices: " << n_verts << " tris: " << n_tris << std::endl;

    // dump obj file
    for (auto& vert : vertices){
        file << "v " << vert.x() << " " << vert.y() << " " << vert.z() << std::endl;
    }
    for (auto& norm : normals){
        file << "vn " << norm.x() << " " << norm.y() << " " << norm.z() << std::endl;
    }
    for (auto& tri : triangles){
        file << "f " << tri[0]+1 << "//" << tri[0]+1 << " " 
                     << tri[1]+1 << "//" << tri[1]+1 << " " 
                     << tri[2]+1 << "//" << tri[2]+1 << std::endl;
    }
    
    GlUtil::Shapes::Sphere* sphere = new GlUtil::Shapes::Sphere();
    sphere->radius = 1;
    auto mesh = GlUtil::DefaultFactories::SphereFactory<GlUtil::MeshTypes::TriangleMesh>::createPrimitive(sphere, 3);
    mesh.toObj("sphere.obj");

    GlUtil::Shapes::TruncatedCone* cone = new GlUtil::Shapes::TruncatedCone();
    cone->radius_bot = 1;
    cone->radius_top = 0.5;
    cone->height = 1;
    cone->fill_top = true;
    cone->fill_bot = true;
    mesh = GlUtil::DefaultFactories::TruncatedConeFactory<GlUtil::MeshTypes::TriangleMesh>::createPrimitive(cone, 10);
    mesh.toObj("cone.obj");
    std::cout << "Cone verts: " << mesh.n_vertices << " cone tris: " << mesh.n_triangles << std::endl;

    // program.addShaderFromSourceFile(QOpenGLShader::Vertex, "C:\\Github\\mCRL2\\tools\\release\\ltsview\\source\\visualizer\\generic_visualizer\\shaders\\simple.vs");
    // program.addShaderFromSourceFile(QOpenGLShader::Fragment, "C:\\Github\\mCRL2\\tools\\release\\ltsview\\source\\visualizer\\generic_visualizer\\shaders\\simple.fs");

    // std::cout << "Program linked: " << (program.link() ? "Succesful. " : "Failed.") << std::endl;
    // std::cout << "Program bound: " << (program.bind() ? "Succesful. " : "Failed.") << std::endl;


    // // Generate vertex arrays and bind
    // glGenVertexArrays(1, &m_VAO);
    // glBindVertexArray(m_VAO);


    
    // total_vertex_data *= 3 * sizeof(float);
    // glGenBuffers(1, &m_PositionBuffer);
    // glBindBuffer(GL_ARRAY_BUFFER, m_PositionBuffer);
    // glBufferData(GL_ARRAY_BUFFER, total_vertex_data, &vertices[0], GL_STATIC_DRAW);
    // glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
    // glEnableVertexAttribArray(0);

    // glGenBuffers(1, &m_NormalBuffer);
    // glBindBuffer(GL_ARRAY_BUFFER, m_NormalBuffer);
    // glBufferData(GL_ARRAY_BUFFER, total_vertex_data, &normals[0], GL_STATIC_DRAW);
    // glVertexAttribPointer(1, 3, GL_FLOAT, GL_TRUE, 0, 0);
    // glEnableVertexAttribArray(1);

    built = true;
}

void GlLTSView::Scene::render(QPainter& painter){
//    m_renderer.render(painter);
}