// Author(s): Rimco Boudewijns and Sjoerd Cranen
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_LTSGRAPH_CAMERA_H
#define MCRL2_LTSGRAPH_CAMERA_H

#include <QQuaternion>
#include <QMatrix4x4>
#include <QVector3D>

/// \brief The arcball is a camera that moves along the surface of
/// three-dimensional sphere of radius 'zoom'. The center
///        of the sphere can move and the camera can move along the surface by
///        applying a rotation. This viewport also contains the width and height
///        of the window. From this it computes the necessary projection and
///        view matrices to draw objects as if they are viewed from this camera.
///        The assumption is that the camera is always upright towards the
///        z-axis. It also uses a perspective projection where objects further
///        away appear smaller.
class ArcballCameraView
{
  public:
  /// \brief Update the view and projection matrix.
  void update();

  /// \brief Converts the given position from world coordinates to Qt window
  /// coordinates.
  QVector3D worldToWindow(QVector3D center) const;

  /// \brief Converts the given position from Qt window coordinates to world
  /// space coordinates.
  QVector3D windowToWorld(QVector3D eye) const;

  /// \brief Set the distance to the camera position.
  void zoom(float zoom)
  {
    m_zoom = zoom;
  }
  float zoom() const
  {
    return m_zoom;
  }

  /// \brief Set the rotation around the camera position.
  void rotation(const QQuaternion& rotation)
  {
    m_rotation = rotation;
  }
  const QQuaternion& rotation() const
  {
    return m_rotation;
  }

  /// \brief Set the position of the camera.
  void center(const QVector3D& position)
  {
    m_center = position;
  }
  const QVector3D& center() const
  {
    return m_center;
  }

  /// \brief Set the dimensions of the viewport.
  void viewport(std::size_t width, std::size_t height)
  {
    m_viewport = QRect(0, 0, static_cast<int>(width), static_cast<int>(height));
  }
  const QRect& viewport() const
  {
    return m_viewport;
  };

  /// \returns The actual position in world coordinates of this camera.
  QVector3D position() const;

  const QMatrix4x4& viewMatrix() const
  {
    return m_viewMatrix;
  }
  const QMatrix4x4& projectionMatrix() const
  {
    return m_projectionMatrix;
  }

  /// \brief Sets the position and orientation of the camera to its default
  /// values.
  void reset();
  /// \brief Sets the rotation of the camera to its default value.
  void resetRotation();

  /// \returns The maximum view distance of this camera.
  float viewdistance() const
  {
    return m_viewdistance;
  }
  void viewdistance(float viewdistance)
  {
    m_viewdistance = viewdistance;
  }

  void translate(const QVector3D& v)
  {
    m_center += v;
  }

  private:
  /// \brief Center of the arcball.
  QVector3D m_center = QVector3D(0.0f, 0.0f, 0.0f);

  /// \brief Rotation around the center.
  QQuaternion m_rotation = QQuaternion(1.0f, 0.0f, 0.0f, 0.0f);

  public:
  /// \brief Distance to the center.
  float m_zoom = 700.0f;

  private:
  float m_viewdistance = 100000.0f;
  float m_vert_fieldofview = 70.0f;

  /// The resulting matrices and viewports, are updated when calling update().
  QMatrix4x4 m_projectionMatrix;
  QMatrix4x4 m_viewMatrix;
  QVector3D m_worldPosition;
  QRect m_viewport;
};

#endif // MCRL2_LTSGRAPH_CAMERA_H
