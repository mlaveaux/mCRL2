// Author(s): Bas Ploeger and Carst Tankink
// Copyright: see the accompanying file COPYING or copy at
// https://github.com/mCRL2org/mCRL2/blob/master/COPYING
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef MATHUTILS_H
#define MATHUTILS_H

namespace MathUtils
{

const double PI = 3.14159265359;

static inline int round_to_int(float f)
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

static inline float deg_to_rad(float deg)
{
  return deg * (float)PI / 180.0f;
}

static inline float rad_to_deg(float rad)
{
  return rad * 180.0f / (float)PI;
}

}
#endif
