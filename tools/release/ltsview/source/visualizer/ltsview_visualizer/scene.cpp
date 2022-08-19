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

#define VERTEX_BINDING 0
#define NORMAL_BINDING 1
#define COLOR_BINDING 2
#define MATRIX_BINDING 3

#define assertGL(msg) msgAssert(reportOpenGLError(), msg)
#define execAssertGL(expr, msg) _result = expr; msgAssert(_result, msg);

bool _result;
inline int cone_encode_resolution(int res, bool a, bool b){
    return res | (a ? (1<<30) : 0) | (b ? (1<<31) : 0);
}
int GlLTSView::ConeDB::getCone(GlUtil::Shapes::TruncatedCone* cone, int resolution, SceneData& sceneData){
    key_t key = {cone->radius_top, cone->radius_bot, cone->height, cone_encode_resolution(resolution, cone->fill_top, cone->fill_bot)};
    auto it = cones.find(key);
    if (it != cones.end()){
        return it->second;
    }else{
        // create the mesh
        GlLTSView::Mesh mesh = GlUtil::DefaultFactories::TruncatedConeFactory<Mesh>::createPrimitive(cone, resolution);
        int mesh_id = (int)sceneData.meshes.size();
        sceneData.meshes.push_back(mesh);
        cones[key] = mesh_id;
        return mesh_id;
    }
}

int GlLTSView::SphereDB::getSphere(GlUtil::Shapes::Sphere* sphere, int resolution, SceneData& sceneData){
    auto it = spheres.find({sphere->radius, resolution});
    if (it != spheres.end()){
        return it->second;
    }else{
        // create the mesh
        GlLTSView::Mesh mesh = GlUtil::DefaultFactories::SphereFactory<Mesh>::createPrimitive(sphere, resolution);
        int mesh_id = (int)sceneData.meshes.size();
        sceneData.meshes.push_back(mesh);
        spheres[{sphere->radius, resolution}] = mesh_id;
        return mesh_id;
    }
}

GlLTSView::Scene::Scene(QOpenGLWidget& glwidget, SceneGraph<GlLTSView::NodeData, GlLTSView::SceneData>& scenegraph, Camera& camera, Cluster *root) : LTSScene(glwidget, scenegraph, camera), m_clusterRoot(root){
}

/// TODO: Fix hardcoded resolutions
template < int sphereRes, int coneRes >
GlLTSView::SGNode* GlLTSView::SceneGraphFunctor<sphereRes, coneRes>::operator()(SGNode* parent, VisTree::VisTreeNode* vistreenode){
    SGNode* node = new SGNode();
    Mesh mesh;
    node->data.data_id = sceneData.matrices.size();
    sceneData.matrices.push_back(vistreenode->data.matrix);
    sceneData.colors.push_back(vistreenode->data.color);
    switch (vistreenode->data.shape->getShapeType()){
        case GlUtil::ShapeType::SPHERE:
            node->data.model_id = sphereDB.getSphere(static_cast<GlUtil::Shapes::Sphere*>(vistreenode->data.shape), sphereRes, sceneData);
            break;
        case GlUtil::ShapeType::TRUNCATED_CONE:
            node->data.model_id = coneDB.getCone(static_cast<GlUtil::Shapes::TruncatedCone*>(vistreenode->data.shape), std::max(coneRes, vistreenode->num_children), sceneData);
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

void GlLTSView::Scene::initialize(){
    if (!built) rebuild();
    initializeOpenGLFunctions();
    gl_funcs = m_glwidget.context()->versionFunctions<QOpenGLFunctions_4_3_Core>();
    
    // create shaderprogram
    execAssertGL(program.addShaderFromSourceFile(QOpenGLShader::Vertex, "shaders/ltsview.vs"), 
                 "Vertex shader can't be added.");
    execAssertGL(program.addShaderFromSourceFile(QOpenGLShader::Fragment,
                                                 "shaders/ltsview.fs"),
                 "Fragment shader can't be added.");

    execAssertGL(program.link(), "Program linking failed.");
    execAssertGL(program.bind(), "Program binding failed.");
    glUseProgram(program.programId());
    genBuffers();
    assertGL("Gen buffers function");
    

    m_view_loc = program.uniformLocation("u_view");
    m_proj_loc = program.uniformLocation("u_proj");
    m_alpha_loc = program.uniformLocation("u_alpha");
}

bool GlLTSView::Scene::reportOpenGLError(){
    GLenum error = glGetError();
    //std::cout << "Checking OGL error." << std::endl;
    if (error != GL_NO_ERROR){
        std::cout << "OpenGL error occured with code 0x" << std::hex << error << ":\n"
                  << program.log().toStdString() << std::endl;
    }
    return error == GL_NO_ERROR;
}

void GlLTSView::Scene::fillBuffers(){
    int index = 0; // used for writing to correct places in arrays
    // Fill vertex SSBO
    // temporarily copy data to gl array
    GLfloat* vertices = (GLfloat*) std::malloc(m_vertices.size() * 4 * sizeof(GLfloat));
    index = 0;
    for (auto& vert : m_vertices){
        vertices[index++] = vert.x();
        vertices[index++] = vert.y();
        vertices[index++] = vert.z();
        vertices[index++] = 1;
    }
    m_vertex_ssbo.bind();
    m_vertex_ssbo.allocate(vertices, index * sizeof(GLfloat));
    assertGL("Vertex ssbo");
    free(vertices);
    

    // same for normals
    // temporarily copy data to gl array
    GLfloat* normals = (GLfloat*) std::malloc(m_normals.size() * 4 * sizeof(GLfloat));
    index = 0;
    for (auto& norm : m_normals){
        normals[index++] = norm.x();
        normals[index++] = norm.y();
        normals[index++] = norm.z();
        normals[index++] = 0;       // important for transforming normals
    }
    m_normal_ssbo.bind();
    m_normal_ssbo.allocate(normals, sizeof(GLfloat) * index);
    assertGL("Normal ssbo");
    free(normals);

    // NB: All of the following buffers should be hinted as dynamic
    // Generate color SSBO
    // temporarily copy data to gl array
    GLfloat* colors = (GLfloat*) std::malloc(m_scenegraph.sceneData.colors.size() * 3 * sizeof(GLfloat));
    index = 0;
    for (auto& color : m_scenegraph.sceneData.colors){
        colors[index++] = color.redF();
        colors[index++] = color.greenF();
        colors[index++] = color.blueF();
    }
    m_color_ssbo.bind();
    m_color_ssbo.allocate(normals, sizeof(GLfloat) * index);
    assertGL("Color ssbo");
    free(colors);

    // finally matrices
    GLfloat* matrices = (GLfloat*) std::malloc(m_scenegraph.sceneData.matrices.size() * 16 * sizeof(GLfloat));
    GLfloat* loc = matrices;
    for (auto& mat : m_scenegraph.sceneData.matrices){
        mat.copyDataTo(loc);
        loc += 16;
    }
    m_matrix_ssbo.bind();
    m_matrix_ssbo.allocate(matrices, m_scenegraph.sceneData.matrices.size() * 16 * sizeof(GLfloat));
    assertGL("Matrix ssbo");
    free(matrices);


    // Fill VBO
    GLuint* vertexdata = (GLuint*) std::malloc(2 * m_vertexData.size() * sizeof(GLuint));
    index = 0;
    for (auto& data : m_vertexData){
        vertexdata[index++] = data.model_index;
        vertexdata[index++] = data.vertex_index;
    }
    m_vbo.bind();
    m_vbo.allocate(vertexdata, index * sizeof(GLuint));
    assertGL("VBO");
    free(vertexdata);

    // tell openGL we will be sending ivec2
    int vertex_location = program.attributeLocation("in_vertexData");
    program.enableAttributeArray(vertex_location);
    program.setAttributeBuffer(vertex_location, GL_FLOAT, 0, 2, 2 * sizeof(GLuint));
    assertGL("Program.enableAttirbuteArray");


    // Fill IBO
    GLuint* indices = (GLuint*) std::malloc(3 * m_triangles.size() * sizeof(GLuint));
    index = 0;
    for (auto& tri : m_triangles){
        indices[index++] = tri[0];
        indices[index++] = tri[1];
        indices[index++] = tri[2];
    }
    m_ibo.bind();
    m_ibo.allocate(indices, index * sizeof(GLuint));
    free(indices);
    assertGL("IBO");
}

void GlLTSView::Scene::genBuffers(){
    assert(built); // we can't create buffers if there's no data
    if (buffers_exist){
        m_vbo.destroy();
        m_ibo.destroy();
        buffers_exist = false;
    }
    assert(!buffers_exist);

    // make sure we are creating buffers in the right place
    m_glwidget.makeCurrent();

    // Create buffers
    m_vbo = QOpenGLBuffer(QOpenGLBuffer::VertexBuffer);
    m_vbo.create();
    assertGL("VBO create");

    m_ibo = QOpenGLBuffer(QOpenGLBuffer::IndexBuffer);
    m_ibo.create();
    assertGL("IBO create");

    m_vertex_ssbo.create();
    m_vertex_ssbo.bind();
    m_vertex_ssbo.setUsagePattern(QOpenGLBuffer::StaticDraw);
    gl_funcs->glBindBufferBase(GL_SHADER_STORAGE_BUFFER, VERTEX_BINDING, m_vertex_ssbo.bufferId());
    assertGL("Vertex SSBO create");

    m_normal_ssbo.create();
    m_normal_ssbo.bind();
    m_normal_ssbo.setUsagePattern(QOpenGLBuffer::StaticDraw);
    gl_funcs->glBindBufferBase(GL_SHADER_STORAGE_BUFFER, NORMAL_BINDING, m_normal_ssbo.bufferId());
    assertGL("Normal SSBO create");

    m_color_ssbo.create();
    m_color_ssbo.bind();
    m_color_ssbo.setUsagePattern(QOpenGLBuffer::DynamicDraw);
    gl_funcs->glBindBufferBase(GL_SHADER_STORAGE_BUFFER, COLOR_BINDING, m_color_ssbo.bufferId());
    assertGL("Color SSBO create");

    m_matrix_ssbo.create();
    m_matrix_ssbo.bind();
    m_matrix_ssbo.setUsagePattern(QOpenGLBuffer::DynamicDraw);
    gl_funcs->glBindBufferBase(GL_SHADER_STORAGE_BUFFER, MATRIX_BINDING, m_matrix_ssbo.bufferId());
    assertGL("Matrix SSBO create");
    
    buffers_exist = true;
}

void GlLTSView::Scene::rebuild(){
    auto start = std::chrono::high_resolution_clock::now();
    // clear the scene data
    m_scenegraph.sceneData.meshes.clear(); m_scenegraph.sceneData.nodes.clear();
    if (built) free(m_scenegraph.root);
    if (built) free(m_vistreeRoot);
    if (built) m_scenegraph.sceneData = SceneData();
    /// TODO: Change settings.h to have enum of modes instead of bool
    m_vistreeRoot = VisTreeGenerator::generate(Settings::instance().clusterVisStyleTubes.value() ? VisTreeGenerator::Mode::TUBES : VisTreeGenerator::Mode::CONES, m_clusterRoot);
    const int sphereRes = 3;
    const int coneRes = 20;
    SceneGraphFunctor<sphereRes, coneRes> sgf(m_scenegraph.sceneData);

    m_scenegraph.root = GlUtil::fold_tree<SGNode, VisTree::VisTreeNode, GlLTSView::SceneGraphFunctor<sphereRes, coneRes>, VisTree::VisTreeNode::childIterator>(
        nullptr, m_vistreeRoot,
        SGNode::setParent, SGNode::addChild,
        sgf, VisTree::VisTreeNode::getChildBegin, VisTree::VisTreeNode::getChildEnd);

    int n_verts = 0;
    int n_tris  = 0; 
    for (const NodeData& node : m_scenegraph.sceneData.nodes){
        n_verts += m_scenegraph.sceneData.meshes[node.model_id].n_vertices;
        n_tris  += m_scenegraph.sceneData.meshes[node.model_id].n_triangles;
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
    for (auto& mesh : m_scenegraph.sceneData.meshes){
        compr_verts += mesh.n_vertices;
        compr_tris += mesh.n_triangles;
    }
    std::cout << "Stored verts: " << compr_verts << ". Stored tris: " << compr_tris << std::endl;
    std::cout << "Total verts: " << n_verts << ". Total tris: " << n_tris << std::endl;

    int vert_offset = 0;
    int tri_offset = 0;
    std::vector<QColor> colors = std::vector<QColor>(n_verts);
    for(auto& mesh : m_scenegraph.sceneData.meshes){
        assert(tri_offset < n_tris);
        assert(vert_offset < n_verts);
    
        std::move(mesh.vertices.begin(), mesh.vertices.end(), back_inserter_vertices);
        std::move(mesh.vertex_normals.begin(), mesh.vertex_normals.end(), back_inserter_normals);
        std::move(mesh.triangles.begin(), mesh.triangles.end(), back_inserter_triangles);

        for (int i = tri_offset; i < tri_offset + mesh.n_triangles; ++i){
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

    std::cout << "Took " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms to complete." << std::endl;
    long uncompr_bytes = n_verts * 96;
    std::cout << "Used resolutions: " << std::endl;
    std::cout << "\t- Sphere: " << sphereRes << std::endl;
    std::cout << "\t- Cone: " << coneRes << " (adaptive)" << std::endl;
    std::cout << uncompr_bytes << " bytes in total for vertices only with no compression with pos/color/transform: " << std::endl;
    std::cout << "\t- " << uncompr_bytes/1024.0 << "KB" << std::endl;
    std::cout << "\t- " << uncompr_bytes/(1024.0*1024.0) << "MB" << std::endl;
    std::cout << "\t- " << uncompr_bytes/(1024.0*1024.0*1024) << "GB" << std::endl;
    long compr_bytes = compr_verts * 16 + m_scenegraph.sceneData.meshes.size() * (64 + 16);
    std::cout << compr_bytes << " bytes in total with compressed storage: " << std::endl;
    std::cout << "\t- " << compr_bytes/1024.0 << "KB" << std::endl;
    std::cout << "\t- " << compr_bytes/(1024.0*1024.0) << "MB" << std::endl;
    std::cout << "\t- " << compr_bytes/(1024.0*1024.0*1024) << "GB" << std::endl;

    auto sphereMesh = GlUtil::DefaultFactories::SphereFactory<GlUtil::Meshes::MeshTypes::TriangleMesh>::createPrimitive(new GlUtil::Shapes::Sphere(), sphereRes);
    sphereMesh.toObj("sphere-res=" + std::to_string(sphereRes) + ".obj");
    auto cone = new GlUtil::Shapes::TruncatedCone();
    cone->radius_top = 1;
    cone->radius_bot = 2;
    cone->height = 1;
    cone->fill_bot = true;
    cone->fill_top = true;
    auto coneMesh = GlUtil::DefaultFactories::TruncatedConeFactory<GlUtil::Meshes::MeshTypes::TriangleMesh>::createPrimitive(cone, coneRes);
    coneMesh.toObj("cone-res=" + std::to_string(coneRes) + ".obj");

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
    // QPainter _painter = QPainter(&m_glwidget);
    std::cout << "Render called" << std::endl;
    painter.beginNativePainting();
    glUseProgram(program.programId());

    /// TODO: Check what needs to be updated.
    QMatrix4x4 view_matrix;
    view_matrix.setToIdentity();
    view_matrix.translate(-10, 0, 0);
    QMatrix4x4 proj_matrix;
    proj_matrix.perspective(45, 16.0/9.0, 0.1, 1000);
    std::cout << "Projection matrix:\n\t";
    for (int i = 0; i < 16; i++){
        std::cout << *(proj_matrix.data() + i);
        if (i%4 == 3) std::cout << std::endl << "\t";
    }
    std::cout << "View matrix:\n\t";
    for (int i = 0; i < 16; i++){
        std::cout << *(view_matrix.data() + i);
        if (i%4 == 3) std::cout << std::endl << "\t";
    }
    // glUniformMatrix4fv(m_view_loc, 1, false, (GLfloat*)m_camera.getViewMatrix().data());
    glUniformMatrix4fv(m_view_loc, 1, false, (GLfloat*)view_matrix.data());
    assertGL("Set uniform; u_view");
    // glUniformMatrix4fv(m_proj_loc, 1, false, (GLfloat*)m_camera.getProjectionMatrix().data());
    glUniformMatrix4fv(m_proj_loc, 1, false, (GLfloat*)proj_matrix.data());
    assertGL("Set uniform; u_proj");

    // set transparency
    glUniform1f(m_alpha_loc, Settings::instance().transparency.value()/100.0f);

    glClearColor(
        Settings::instance().backgroundColor.value().red() / 255.0,
        Settings::instance().backgroundColor.value().green() / 255.0,
        Settings::instance().backgroundColor.value().blue() / 255.0,
        1.0f
    );
    // Clear existing buffers
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS); // closer to camera wins

    glCullFace(GL_FRONT);
    assertGL("GlCullFace(GL_FRONT)");


    /// TODO: Multiple calls to get proper transparency.
    std::cout << "Before draw elements" << std::endl;
    m_vbo.bind();
    m_ibo.bind();
    glDrawElements(GL_TRIANGLES, m_triangles.size() * 3, GL_UNSIGNED_INT, nullptr);
    std::cout << "After draw elements" << std::endl;
    assertGL("GlDrawElements for back faces");

    glClear(GL_DEPTH_BUFFER_BIT);
    glCullFace(GL_BACK);
    assertGL("GlCullFace(GL_BACK)");

    /// TODO: Multiple calls to get proper transparency.
    glDrawElements(GL_TRIANGLES, m_triangles.size() * 3, GL_UNSIGNED_INT, nullptr);
    assertGL("GlDrawElements for front faces");

    painter.endNativePainting();
}