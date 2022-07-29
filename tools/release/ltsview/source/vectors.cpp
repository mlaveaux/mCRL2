// Author(s): Bas Ploeger
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
/// \file vectors.cpp
/// \brief Implementation of the Vector classes

#include <cmath>
#include "mathutils.h"
#include "vectors.h"

namespace Vectors {

void toPolar(float &angle, float &radius, const QVector2D& v)
{
  angle = MathUtils::rad_to_deg(atan2(v.y(), v.x()));
  if (angle < 0.0f)
  {
    angle += 360.0f;
  }
  radius = v.length();
}

QVector2D fromPolar(float angle, float radius)
{
  float angle_rad = MathUtils::deg_to_rad(angle);
  return QVector2D(radius * cos(angle_rad), radius * sin(angle_rad));
}

}
