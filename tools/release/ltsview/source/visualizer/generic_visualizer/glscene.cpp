// Author(s): Ruben Vink.
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
// Abstracted from ltsgraph/glscene.cpp (02/08/2022)

#include "glscene.h"
#include "glcamera.h"
#include "glrenderer.h"
#include "glscenegraph.h"

template<typename NodeData, typename ModelData>
void GLScene<NodeData, ModelData>::initialize(){
    initializeOpenGLFunctions();

    m_scenegraph.rebuild();
    m_renderer.initialize(m_scenegraph);
}

/// \brief Updates the state of the scene.
template<typename NodeData, typename ModelData>
void GLScene<NodeData, ModelData>::update(){
    m_camera.update();
    m_scenegraph.update(); // often the scenegraph and renderer do not actually
    m_renderer.update();   // require update since the image is static
}

template<typename NodeData, typename ModelData>
void GLScene<NodeData, ModelData>::rebuild(){
    m_camera.reset();
    m_scenegraph.rebuild();
    m_renderer.rebuild(m_scenegraph);
}

template<typename NodeData, typename ModelData>
void GLScene<NodeData, ModelData>::resize(std::size_t width, std::size_t height){
    if (m_camera){
        m_camera->setViewport(width, height);
    }
}

template<typename NodeData, typename ModelData>
void GLScene<NodeData, ModelData>::render(QPainter& painter){
    // we let the renderer do it's magic
    m_renderer->render(painter);
}

