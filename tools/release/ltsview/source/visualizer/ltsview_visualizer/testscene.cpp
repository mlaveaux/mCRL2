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


const char* VS = "#version 430 compatibility\n"

"layout(std140) buffer b_test { vec4 colors[]; };\n"

"layout(location = 0) in ivec3 vert;\n"
"out vec4 vColor;\n"


"uniform int u_n_colors;\n"

"uniform int u_offset_colors;\n"

"void main(void)\n"
"{\n"
"   gl_Position = vec4(vert.xy/10.0, 0, 1);\n"
"   vColor      = vec4(1);\n"
"}";

const char* FS = "#version 430 compatibility\n"

"in vec4 vColor;\n"
"out vec4 fColor;\n"
"void main(void)\n"
"{\n"
"   fColor = vColor;\n"
"}";

Test::TScene::TScene(Cluster* root) : TestScene(){
    m_clusterRoot = root;
}

/// @brief Initialize scene. Assumption that initializeOpenGLFunctions() has been called in advance.
void Test::TScene::initializeScene(){
    m_prog = new QOpenGLShaderProgram;

    m_prog->addShaderFromSourceCode(QOpenGLShader::Vertex, VS);
    m_prog->addShaderFromSourceCode(QOpenGLShader::Fragment, FS);
    m_prog->bindAttributeLocation("vertex", 0);
    m_prog->bindAttributeLocation("fragment", 1);
    m_prog->link();

    m_prog->bind();
    
    
    mCRL2log(mcrl2::log::debug) << "Program linked and bound." << std::endl;

    QOpenGLFunctions *f = QOpenGLContext::currentContext()->functions();

    m_vao.create();
    QOpenGLVertexArrayObject::Binder vaoBinder(&m_vao);

    //VBO
    m_vbo.create();
    m_vbo.bind();
        m_vbo.allocate(new int[9]{-10, -10, 0, 10, -10, 0, 0, 10, 0}, 9*sizeof(int));
        f->glEnableVertexAttribArray(0);
        f->glVertexAttribPointer(0, 3, GL_INT, GL_FALSE, 0, 0);
    m_vbo.release();

    //IBO
    // ogl->glGenBuffers(1, &m_ibo);
    // ogl->glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_ibo);
    // ogl->glBufferData(GL_ELEMENT_ARRAY_BUFFER, 3*4, new int[3]{0, 2, 1}, GL_STATIC_DRAW);

    // //SSBO
    // glGenBuffers(1, &m_ssbo);
    // // Fill buffer
    // glBindBuffer(GL_SHADER_STORAGE_BUFFER, m_ssbo);
    // glBufferData(GL_SHADER_STORAGE_BUFFER, 3*4*4, new float[12]{1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1}, GL_STATIC_DRAW);
    // glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT); /// TODO: test if necessary
    // glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
    // // Tell shaders where to find this buffer
    // GLuint loc = glGetProgramResourceIndex(program.programId(),
    //                                         GL_SHADER_STORAGE_BLOCK, "b_test");
    // // Bind buffer base
    // glBindBufferBase(GL_SHADER_STORAGE_BUFFER, loc, m_ssbo);

    u_n_colors = m_prog->uniformLocation("u_n_colors");
    u_offset_colors = m_prog->uniformLocation("u_offset_colors");
    glCheckError();

    mCRL2log(mcrl2::log::debug) << "[testscene] Device name: " << f->glGetString(GL_RENDERER) << std::endl;
}


void Test::TScene::renderScene()
{
    QOpenGLExtraFunctions *f = QOpenGLContext::currentContext()->extraFunctions();
    f->glClear(GL_COLOR_BUFFER_BIT);
    
    int* xy = new int[4];
    f->glGetIntegerv(GL_VIEWPORT, xy);
    mCRL2log(mcrl2::log::debug) << "renderScene viewport sizes: (" << xy[0] << ", " <<
                                                                      xy[1] << ", " <<
                                                                      xy[2] << ", " <<
                                                                      xy[3] << ")\n"; 

    m_prog->bind();

    QOpenGLVertexArrayObject::Binder vaoBinder(&m_vao);
    int value;
    f->glGetIntegerv(GL_CONTEXT_FLAGS, &value);
    // mCRL2log(mcrl2::log::debug) << "GL_CONTEXT_FLAGS: " << QString("%1").arg(value, 16, 2, QChar('0')).toStdString() << std::endl; 
    int major, minor;
    f->glGetIntegerv(GL_MAJOR_VERSION, &major);
    f->glGetIntegerv(GL_MINOR_VERSION, &minor);
    mCRL2log(mcrl2::log::debug) << "GL Version: " << major << "." << minor << std::endl;
    // mCRL2log(mcrl2::log::debug) << "GL_CONTEXT_FLAG_FORWARD_COMPATIBLE_BIT: " << GL_CONTEXT_FLAG_FORWARD_COMPATIBLE_BIT << std::endl;
    // mCRL2log(mcrl2::log::debug) << "GL_CONTEXT_FLAG_DEBUG_BIT: " << GL_CONTEXT_FLAG_DEBUG_BIT << std::endl;
    // mCRL2log(mcrl2::log::debug) << "GL_CONTEXT_FLAG_ROBUST_ACCESS_BIT: " << GL_CONTEXT_FLAG_ROBUST_ACCESS_BIT << std::endl;
    // mCRL2log(mcrl2::log::debug) << "GL_CONTEXT_FLAG_NO_ERROR_BIT: " << GL_CONTEXT_FLAG_NO_ERROR_BIT << std::endl;
    int num_versions;


    // f->glGetIntegerv(GL_NUM_SHADING_LANGUAGE_VERSIONS, &num_versions);
    // mCRL2log(mcrl2::log::debug) << "Found " << num_versions << " supported shading language versions." << std::endl;
    // const char* version_name = new char[128];
    // for (int k = 0; k < num_versions; k++){
    //     mCRL2log(mcrl2::log::debug) << '\t' << f->glGetStringi(GL_SHADING_LANGUAGE_VERSION, k) << std::endl;
    // }

    f->glDrawArrays(GL_TRIANGLES, 0, 9);
    glCheckError();
    vaoBinder.release();
    m_prog->release();
}