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

"layout(std430) buffer b_test { vec4 colors[]; };\n"

"layout(std430) buffer b_alpha { float alpha[]; };\n"

"layout(location = 0) in ivec3 vert;\n"
"out vec4 vColor;\n"


"uniform int u_n_colors;\n"

"uniform int u_offset_colors;\n"

"void main(void)\n"
"{\n"
"   gl_Position = vec4(vert.xy/10.0, 0, 1);\n"
"   int offset  = (u_offset_colors + vert.z) % u_n_colors;\n"
"   vColor      = vec4(vec3(colors[offset]), alpha[offset]);\n"
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
    m_prog->link();

    m_prog->bind();
    
    
    mCRL2log(mcrl2::log::debug) << "Program linked and bound." << std::endl;

    QOpenGLExtraFunctions *f = QOpenGLContext::currentContext()->extraFunctions();
    QOpenGLFunctions_4_3_Core *f430 = QOpenGLContext::currentContext()->versionFunctions<QOpenGLFunctions_4_3_Core>();
    m_vao.create();
    QOpenGLVertexArrayObject::Binder vaoBinder(&m_vao);

    //VBO
    m_vbo.create();
    m_vbo.bind();
        m_vbo.allocate(new int[9]{-10, -10, 0, 10, -10, 1, 0, 10, 2}, 9*sizeof(int));
        f->glEnableVertexAttribArray(0);
        f->glVertexAttribIPointer(0, 3, GL_INT, 0, 0);
    m_vbo.release();

    //IBO
    m_ibo = QOpenGLBuffer(QOpenGLBuffer::IndexBuffer);
    m_ibo.create();
    m_ibo.bind();
        m_ibo.allocate(new int[3]{0, 2, 1}, 3*sizeof(int));
    m_ibo.release();

    //SSBO
    int ssbo_count = 0;
    m_ssbo1.create();
    m_ssbo1.bind();
        float alpha = 0;
        m_ssbo1.allocate(new float[12]{1, 0, 0, alpha, 0, 1, 0, alpha, 0, 0, 1, alpha}, 12*sizeof(float));
        GLuint SSBO1_block_index = f->glGetProgramResourceIndex(m_prog->programId(), GL_SHADER_STORAGE_BLOCK, "b_test");
        f430->glShaderStorageBlockBinding(m_prog->programId(), SSBO1_block_index, ssbo_count);
        f->glBindBufferBase(GL_SHADER_STORAGE_BUFFER, ssbo_count++, m_ssbo1.bufferId());
    m_ssbo1.release();

    m_ssbo2.create();
    m_ssbo2.bind();
        m_ssbo2.allocate(new float[3]{1, 0.6, 0.2}, 3*sizeof(float));
        GLuint SSBO2_block_index = f->glGetProgramResourceIndex(m_prog->programId(), GL_SHADER_STORAGE_BLOCK, "b_alpha");
        f430->glShaderStorageBlockBinding(m_prog->programId(), SSBO2_block_index, ssbo_count);
        f->glBindBufferBase(GL_SHADER_STORAGE_BUFFER, ssbo_count++, m_ssbo2.bufferId());
    m_ssbo2.release();

    m_prog->setUniformValue("u_n_colors", 3);
    glCheckError();

    mCRL2log(mcrl2::log::debug) << "[testscene] Device name: " << f->glGetString(GL_RENDERER) << std::endl;
}

int color_offset = 0;
void Test::TScene::renderScene()
{
    QOpenGLExtraFunctions *f = QOpenGLContext::currentContext()->extraFunctions();
    f->glEnable(GL_BLEND);
    f->glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    int* xy = new int[4];
    f->glGetIntegerv(GL_VIEWPORT, xy);
    mCRL2log(mcrl2::log::debug) << "renderScene viewport sizes: (" << xy[0] << ", " <<
                                                                      xy[1] << ", " <<
                                                                      xy[2] << ", " <<
                                                                      xy[3] << ")\n"; 

    m_prog->bind();
    m_prog->setUniformValue("u_offset_colors", color_offset++);

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
    m_ibo.bind();
    // m_ssbo1.bind();
    f->glDrawElements(GL_TRIANGLES, 3, GL_UNSIGNED_INT, (void*)0);
    glCheckError();
    vaoBinder.release();
    m_prog->release();
}