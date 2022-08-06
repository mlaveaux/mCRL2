// Author(s): Ruben Vink.
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
// Abstracted from ltsgraph/glscene.h (02/08/2022)

#ifndef MCRL2_OPENGL_GLSCENE_H
#define MCRL2_OPENGL_GLSCENE_H

#include <array>
#include <QOpenGLBuffer>
#include <QOpenGLFunctions_3_3_Core>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLWidget>
#include <QPainter>

#include "glscenegraph.h"
#include "glrenderer.h"

class Camera;

template<typename NodeData, typename ModelData>
class GLScene : private QOpenGLFunctions_3_3_Core
{
public:
  /// \brief Constructor.
  /// \param glwidget The widget where this scene is drawn
  /// \param scenegraph The scenegraph that is to be drawn.
  /// \param cam Camera to be used
  GLScene(QOpenGLWidget& glwidget, SceneGraph<NodeData, ModelData>& scenegraph, Camera& camera, Renderer& renderer) :
   m_glwidget(glwidget),
   m_scenegraph(scenegraph),
   m_camera(camera),
   m_renderer(renderer) {}

  /// \brief Initializes all state and data required for rendering.
  void initialize();

  /// \brief Updates the state of the scene.
  void update();

  /// \brief Builds static data based on the current data.
  void rebuild();

  /// \brief Render the scene.
  void render(QPainter& painter);

  /// \brief Resize the OpenGL viewport.
  void resize(std::size_t width, std::size_t height);

private:
  Camera& m_camera;
  SceneGraph<NodeData, ModelData>& m_scenegraph;
  Renderer& m_renderer;
  QOpenGLWidget& m_glwidget; /// The widget where this scene is drawn

  float m_device_pixel_ratio;
  QFont m_font;
};
#endif // MCRL2_OPENGL_GLSCENE_H
