#include "glluts.h"
#include <cmath>

using namespace GlUtil::LUTs;

std::vector<float> CircleLUT::LUTcosf = {};
std::vector<float> CircleLUT::LUTsinf = {};
int CircleLUT::last_resolution = -1;


void CircleLUT::update(int resolution)
{
  if (last_resolution != resolution)
    recompute(resolution);
}

void CircleLUT::recompute(int resolution)
{
  
  LUTcosf = std::vector<float>(resolution); // allocate space for arrays
  LUTsinf = std::vector<float>(resolution); // allocate space for arrays
  
  float a = 0, da = 6.28318530718 / resolution;
  for (int i = 0; i < resolution; ++i)
  {
    LUTcosf[i] = std::cosf(a);
    LUTsinf[i] = std::sinf(a);
    a += da;
  }
}