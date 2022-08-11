#ifndef MCRL2_OPENGL_LUTS_H
#define MCRL2_OPENGL_LUTS_H

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
} // namespace LUTs
} // namespace GlUtil
#endif