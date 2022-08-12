// Author(s): Ruben Vink.
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_OPENGL_CAMERA_H
#define MCRL2_OPENGL_CAMERA_H

#include <QVector3D>
#include <QMatrix4x4>

/// \brief Camera is an abstract class defining what cameras should implement
///        in order to work with the OpenGL structure.
class Camera
{
  public:
  /// \brief Converts the given position from world space coordinates to
  /// projection space coordinates. Uses OpenGL coordinates \param ws_vector
  /// WorldSpace_vector \returns Position of ws_vector in projection space
  QVector3D worldToWindow(QVector3D ws_vector) const
  {
    return ws_vector.project(m_viewMatrix, m_projectionMatrix, m_viewport);
  }

  /// \brief Converts the given position from projection space coordinates to
  /// world space coordinates.
  QVector3D windowToWorld(QVector3D proj_vector) const
  {
    return proj_vector.unproject(m_viewMatrix, m_projectionMatrix, m_viewport);
  }

  // Getters for various properties
  QVector3D getPosition() const { return m_worldPosition; }
  const QMatrix4x4& getViewMatrix() const { return m_viewMatrix; }
  const QMatrix4x4& getProjectionMatrix() const { return m_projectionMatrix; }
  const QRect& getViewport() const { return m_viewport; }

  // setters for various properties
  void setPosition(const QVector3D& position) { m_worldPosition = position; }
  void setViewport(std::size_t width, std::size_t height) { m_viewport = QRect(0, 0, (int)width, (int)height); }

  // virtual functions
  virtual void reset() = 0;
  virtual void update() = 0;
  protected:
  // Matrices to be sent to GPU for rendering
  QMatrix4x4 m_projectionMatrix; // view  -> projection
  QMatrix4x4 m_viewMatrix;       // world -> view
  // Properties of camera
  QVector3D m_worldPosition;
  QRect m_viewport;
};

#endif