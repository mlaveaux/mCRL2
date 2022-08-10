#ifndef MCRL2_LTSVIEW_RENDERER_H
#define MCRL2_LTSVIEW_RENDERER_H

#include "scene.h"
#include "glutil.h"

class LTSRenderer : private QOpenGLFunctions_3_3_Core{
    public:
        void render(QPainter& painter);

        void initializeShaders();
        void initializeData(SceneGraph<GlLTSView::NodeData, GlLTSView::SceneData>& sg);

    private:
        GLuint m_VAO;
        GLuint m_PositionBuffer;
        GLuint m_NormalBuffer;
        QOpenGLShaderProgram program;
};

#endif