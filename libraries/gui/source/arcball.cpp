// Author(s): Johannes Altmanninger
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include "mcrl2/gui/arcball.h"

#include "mcrl2/gui/glu.h"

#ifdef MCRL2_PLATFORM_MAC
  #include <GLKit/GLKMatrix4.h>
#endif

#ifndef M_PI
  #define M_PI 3.14159265358979323846264338327950288
#endif

namespace mcrl2
{
namespace gui
{

/// \brief Compute the projection of the point p on a sphere that covers the viewport
static QVector3D arcballVector(const QPoint& p)
{
  float x, y, z;
  const float radius = 1.5f;
  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT, viewport);

  // Compute x and y relative to the middle of the viewport
  x = (float)p.x() / viewport[2] * 2.0f - 1.0f;
  y = (float)p.y() / viewport[3] * 2.0f - 1.0f;
  y = -y;

  float squared = x * x + y * y;
  if (squared <= radius * radius)
  {
    z = std::sqrt(radius * radius - squared);
  }
  else
  {
    float len = std::sqrt(squared);
    if (std::isnormal(len))
    {
      x /= len;
      y /= len;
    }
    else
    {
      x = 0.0f;
      y = 0.0f;
    }
    z = 0.0f;
  }
  return QVector3D(x, y, z).normalized();
}

QQuaternion arcballRotation(const QPoint& p1, const QPoint& p2, const float rotation_speed)
{
  // update rotation based on the difference in mouse coordinates
  QVector3D v1 = arcballVector(p1).normalized();
  QVector3D v2 = arcballVector(p1 + rotation_speed * (p2 - p1)).normalized();
  QVector3D cross = QVector3D::crossProduct(v1, v2);
  float dot = QVector3D::dotProduct(v1, v2);
  return QQuaternion(dot, cross);
}

void applyRotation(const QQuaternion& rotation, bool reverse)
{
  float angle = 180 / M_PI * std::acos(std::min(1.0f, rotation.scalar()));
  if (reverse)
  {
    angle = -angle;
  }
  // not sure why the angle has to be doubled
  glRotatef(2 * angle, rotation.x(), rotation.y(), rotation.z());
}

} // namespace gui
} // namespace mcrl2
