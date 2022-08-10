#ifndef MCRL2_OPENGL_GLLUTS_H
#define MCRL2_OPENGL_GLLUTS_H

#include <cstdlib>
#include <vector>

namespace GlUtil
{
namespace LUTs
{
class CircleLUT
{
  public:
  static std::vector<float> LUTcosf;
  static std::vector<float> LUTsinf;
  static int last_resolution;
  static void update(int resolution);

  private:
  CircleLUT()
  {
  }
  static void recompute(int resolution);
};

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
    LUTcosf[i] = cosf(a);
    LUTsinf[i] = sinf(a);
    a += da;
  }
}
} // namespace LUTs
} // namespace GlUtil
#endif