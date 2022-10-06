// Author(s): Rimco Boudewijns and Sjoerd Cranen
// Generalized by Ruben Vink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include "arcballcamera.h"

#include "mcrl2/gui/arcball.h"
#include "mcrl2/gui/glu.h"
#include "mcrl2/utilities/logger.h"

QQuaternion ArcballCamera::DEFAULT_ROTATION = QQuaternion(1.0f, 0.0f, 0.0f, 0.0f);
QVector3D ArcballCamera::DEFAULT_CENTER   = QVector3D(0, 0, 0);

void ArcballCamera::update()
{
  m_viewMatrix.setToIdentity();
  m_viewMatrix.translate(QVector3D(0.0f, 0.0f, -m_zoom));
  m_viewMatrix.rotate(m_rotation);
  m_viewMatrix.translate(m_center);

  QMatrix4x4 invertedViewMatrix = m_viewMatrix.inverted();
  m_worldPosition = QVector3D(invertedViewMatrix.column(3));

  m_projectionMatrix.setToIdentity();
  m_projectionMatrix.perspective(m_vert_fieldofview, static_cast<float>(m_viewport.width()) / static_cast<float>(m_viewport.height()), 1.0f, m_viewdistance);

  mCRL2log(mcrl2::log::debug) << "Cam settings; zoom: " << m_zoom 
  << " worldpos: (" << m_worldPosition.x() << ", " << m_worldPosition.y() << ", " << m_worldPosition.z() << ") center: (" << m_center.x() << ", " << m_center.y() << ", " << m_center.z() << ") FOV: " << 
  m_vert_fieldofview << " viewport: (" << m_viewport.width() << ", " << m_viewport.height() << ")" << std::endl;

}

void ArcballCamera::applyTranslate(const QPoint& p1, const QPoint& p2, const float translation_speed){
  QVector3D translation(windowToWorld(QVector3D(p2.x(), p2.y(), m_worldPosition.z())) - windowToWorld(QVector3D(p1.x(), p1.y(), m_worldPosition.z())));
  m_center += translation_speed*translation;
  update();
}

void ArcballCamera::applyDragRotate(const QPoint& p1, const QPoint& p2, const float rotation_speed){
  m_rotation = mcrl2::gui::arcballRotation(p1, p2) * m_rotation;
  update();
}

void ArcballCamera::applyZoom(const QPoint& p1, const QPoint& p2, const float zoom_speed, const float zoom_base){
  QPoint diff = p2 - p1;
  m_zoom *= std::powf(zoom_base, diff.y()*zoom_speed);
  m_zoom = std::max(m_zoom, 0.0001f); // if zoom gets to 0 we don't have zoom anymore
  update();
}

void ArcballCamera::reset()
{
  m_center = ArcballCamera::DEFAULT_CENTER;
  m_rotation = ArcballCamera::DEFAULT_ROTATION;
  m_zoom = ArcballCamera::DEFAULT_ZOOM;
}

void ArcballCamera::resetRotation()
{
  m_rotation = ArcballCamera::DEFAULT_ROTATION;
}
