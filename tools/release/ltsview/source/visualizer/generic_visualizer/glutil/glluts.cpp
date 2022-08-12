#include "glluts.h"
#include <cmath>
#include "mathutils.h"

using namespace GlUtil::LUTs;

std::vector<qreal> CircleLUT::LUTcosf = {};
std::vector<qreal> CircleLUT::LUTsinf = {};
int CircleLUT::last_resolution = -1;


void CircleLUT::update(int resolution)
{
  if (last_resolution != resolution)
    recompute(resolution);
}

void CircleLUT::recompute(int resolution)
{
  
  LUTcosf = std::vector<qreal>(resolution); // allocate space for arrays
  LUTsinf = std::vector<qreal>(resolution); // allocate space for arrays
  
  qreal a = 0, da = MathUtils::TWO_PI / resolution;
  for (int i = 0; i < resolution; ++i)
  {
    LUTcosf[i] = std::cos(a);
    LUTsinf[i] = std::sin(a);
    a += da;
  }
}