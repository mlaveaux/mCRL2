// Author(s): Rimco Boudewijns and Sjoerd Cranen
// Generalized by Ruben Vink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef MCRL2_OPENGL_ARCBALLCAMERA_H
#define MCRL2_OPENGL_ARCBALLCAMERA_H

#include "glcamera.h"
#include <QQuaternion>
#include <QMatrix4x4>
#include <QVector3D>

/// \brief The arcball is a camera that moves along the surface of three-dimensional sphere of radius 'zoom'. The center
///        of the sphere can move and the camera can move along the surface by applying a rotation. This viewport
///        also contains the width and height of the window. From this it computes the necessary projection
///        and view matrices to draw objects as if they are viewed from this camera. The assumption is that the camera
///        is always upright towards the z-axis. It also uses a perspective projection where objects further away appear
///        smaller.
class ArcballCamera : public Camera
{
public:
  /// \brief Update the view and projection matrix.
  void update() override;

  /// \brief Set the distance to the camera position.
  void setZoom(float zoom) { m_zoom = zoom; }
  float getZoom() const { return m_zoom; }

  /// \brief Set the rotation around the camera position.
  void setRotation(const QQuaternion& rotation) { m_rotation = rotation; }
  const QQuaternion& getRotation() const { return m_rotation; }

  /// \brief Set the position of the camera.
  void setCenter(const QVector3D& center) { m_center = center; }
  const QVector3D& getCenter() const { return m_center; }

  /// \brief Sets the position and orientation of the camera to its default values.
  void reset();
  /// \brief Sets the rotation of the camera to its default value.
  void resetRotation();

  /// \returns The maximum view distance of this camera.
  float getViewdistance() const { return m_farPlane; }
  void setViewdistance(float viewdistance) { m_farPlane = viewdistance; }

  /// \brief Rotate camera around center using drag coordinates.
  /// \param p1 Old point
  /// \param p2 New point
  /// \param rotation_speed TODO: implement user control
  void applyDragRotate(const QPoint& p1, const QPoint& p2, const float rotation_speed = 1.0f);

  /// @brief Applies zoom towards center using drag coordinates.
  /// @param p1 Old point
  /// @param p2 New point
  /// @param zoom_speed How fast should it zoom (linear)
  /// @param zoom_base How fast should it zoom (exponential)
  void applyZoom(const QPoint& p1, const QPoint& p2, const float zoom_speed = 0.5f, const float zoom_base = 1.005);

  /// @brief Applies translation of the center using drag coordinates.
  /// @param p1 Old point
  /// @param p2 New point
  /// @param translation_speed How fast should it translate (linear) 
  void applyTranslate(const QPoint& p1, const QPoint& p2, const float translation_speed = 10.0f);

private:
  /// \brief Center of the arcball.
  QVector3D m_center = ArcballCamera::DEFAULT_CENTER;

  /// \brief Rotation around the center.
  QQuaternion m_rotation = ArcballCamera::DEFAULT_ROTATION;

  /// \brief Distance to the center.
  float m_zoom = ArcballCamera::DEFAULT_ZOOM;

  float m_nearPlane = ArcballCamera::DEFAULT_NEAR_PLANE;
  float m_farPlane = ArcballCamera::DEFAULT_FAR_PLANE;
  float m_vert_fieldofview = ArcballCamera::DEFAULT_VERTFIELDOFVIEW;
public:
  // Default values
  static constexpr float DEFAULT_VERTFIELDOFVIEW = 75.0f;
  static constexpr float DEFAULT_ZOOM = 1.0f;
  static constexpr float DEFAULT_NEAR_PLANE = 0.01f;
  static constexpr float DEFAULT_FAR_PLANE = 10000.0f;
  static QQuaternion DEFAULT_ROTATION;
  static QVector3D DEFAULT_CENTER;
};

#endif // MCRL2_LTSGRAPH_CAMERA_H
