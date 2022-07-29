// Author(s): Ruben Vink.
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_OPENGL_RENDERER_H
#define MCRL2_OPENGL_RENDERER_H

#include "glscenegraph.h"
#include <QOpenGLShaderProgram>
#include <QPainter>

/// \brief Objects of type renderer represent a pipeline of shaders.
template <typename NodeData, typename ModelData>
class Renderer : private QOpenGLShaderProgram
{
    public:
        void initialize(SceneGraph<NodeData, ModelData>* scenegraph);
        void forceInitialize(SceneGraph<NodeData, ModelData>* scenegraph);

        virtual void update();
        // in general usage as follows:
        // painter.beginNativePainting()
        // ~ OpenGL code
        // painter.endNativePainting()
        // ~ Draw some text using QPainter
        virtual void render(QPainter& painter) = 0;
        // used to clear buffers etc.
        virtual void rebuild(SceneGraph<NodeData, ModelData>* scenegraph) = 0;


    private:
        virtual void intializeShaders();
        virtual void initializeData(SceneGraph<NodeData, ModelData>* scenegraph);
        bool m_initialized = false;
};

#endif // MCRL2_OPENGL_RENDERER_H