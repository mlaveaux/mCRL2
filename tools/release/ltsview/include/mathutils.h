// Author(s): Bas Ploeger and Carst Tankink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef MATHUTILS_H
#define MATHUTILS_H
#include <QtCore>

namespace MathUtils
{

const inline qreal PI = 3.14159265359;
const inline qreal TWO_PI = 2*PI;
const inline qreal INV_PI = 1/PI;
const inline qreal INV_SQRT2 = 0.70710678118;
const inline qreal INV_SQRT2F = 0.70710678118; 
const inline qreal ONE_THIRD  = 0.33333333333;

static inline int round_to_int(qreal f)
{
  if (f >= 0.0f)
  {
    return (int)(f + 0.5f);
  }
  else
  {
    return (int)(f - 0.5f);
  }
}

const inline qreal deg2rad = PI/180;
static inline qreal deg_to_rad(qreal deg)
{
  return deg * deg2rad;
}

const inline qreal rad2deg = 180 * INV_PI;
static inline qreal rad_to_deg(qreal rad)
{
  return rad * rad2deg;
}

}
#endif
