#ifndef MCRL2_OPENGL_PRIMITIVEFACTORIES_H
#define MCRL2_OPENGL_PRIMITIVEFACTORIES_H

#include "glshapes.h"
#include "glmeshes.h"
#include <cstdlib>
#include <cstring>
#include "glluts.h"

#include <map>
#include <algorithm>

namespace GlUtil
{
// The goal is to have a primitive factory template class that can be easily
// wrapped inside a functor to
//   map over a VisTree
template <typename T, typename Shape, typename MeshType>
class AbstractPrimitiveFactory
{
  public:
  static MeshType createPrimitive(Shape* shape, int resolution);

  private:
  AbstractPrimitiveFactory()
  {
  }
};

namespace DefaultFactories
{

template <typename MeshType>
class RingFactory : AbstractPrimitiveFactory<RingFactory<MeshType>,
                                               Shapes::Ring, MeshType>
{
  public:
  static MeshType createPrimitive(Shapes::Ring* shape,
                                  int resolution);
};


template <typename MeshType>
class SphereFactory : AbstractPrimitiveFactory<SphereFactory<MeshType>,
                                               Shapes::Sphere, MeshType>
{
  public:
  static MeshType createPrimitive(Shapes::Sphere* shape,
                                  int resolution);
};

template <typename MeshType>
class TruncatedConeFactory
    : AbstractPrimitiveFactory<TruncatedConeFactory<MeshType>,
                               Shapes::TruncatedCone, MeshType>
{
  public:
  static MeshType createPrimitive(Shapes::TruncatedCone* shape,
                                  int resolution);

  private:
  static MeshTypes::Vertices last_ring;
  static int last_resolution;
};

} // namespace DefaultFactories
} // namespace GlUtil
#endif