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
#include <QOpenGLWidget>
#include <QPainter>

#include "glscenegraph.h"
#include "glcamera.h"

template<typename NodeData, typename SceneData>
class GLScene
{
public:
  GLScene(){};

  /// \brief Initializes all state and data required for rendering.
  virtual void initializeScene() = 0;

  /// \brief Updates the state of the scene.
  virtual void updateScene() = 0;

  /// \brief Builds static data based on the current data.
  virtual void rebuildScene() = 0;

  /// \brief Render the scene using pure OpenGL.
  virtual void renderScene() = 0;

  /// \brief Resize the OpenGL viewport.
  virtual void resizeScene(std::size_t width, std::size_t height) = 0;

protected:
  Camera* m_camera;
  SceneGraph<NodeData, SceneData> m_scenegraph;
  float m_device_pixel_ratio;
  QFont m_font;
};
#endif // MCRL2_OPENGL_GLSCENE_H
