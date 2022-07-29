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
