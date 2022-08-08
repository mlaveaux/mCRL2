#ifndef MCRL2_OPENGL_PRIMITIVES_H
#define MCRL2_OPENGL_PRIMITIVES_H

#include <QMatrix4x4>

// The goal is to have a primitive factory template class that can be easily
// wrapped inside a functor to
//   map over a VisTree
template <typename T, typename Shape, typename ModelData>
class AbstractPrimitiveFactory
{
  public:
  static ModelData createPrimitive(Shape& shape, int resolution);

  private:
  AbstractPrimitiveFactory()
  {
  }
};

namespace Primitives
{

enum ShapeType
{
  SPHERE,
  HEMISPHERE,
  TUBE,
  CONE,
  TRUNCATED_CONE,
  OBLIQUE_CONE,
  DISC,
  RING,
};

struct TriangleMesh
{
  int n_vertices;
  float* vertices;
  float* vertex_normals;
  int n_triangles;
  int* triangles;
  float* face_normals;
};

struct QuadMesh
{
  int n_vertices;
  float* vertices;
  float* vertex_normals;
  int n_quads;
  int* quads;
  float* face_normals;
};

struct Shape
{
  ShapeType shape;
};

namespace Shapes
{
struct Sphere : Shape
{
  ShapeType shape = ShapeType::SPHERE;
  float radius;
};

struct HemiSphere : Shape
{
  ShapeType shape = ShapeType::HEMISPHERE;
  float radius;
};

// without transform a tube is formed by circle in the xy0 and xy1 planes
struct Tube : Shape
{
  ShapeType shape = ShapeType::TUBE;
  float radius;
  float height;
};

// without transform a cone is formed by a cirle in the xy0 plane and a point at
// 001
struct Cone : Shape
{
  ShapeType shape = ShapeType::CONE;
  float radius;
  float height;
};

// without transform a disc is formed by a circle in the xy0 plane (filled)
struct Disc : Shape
{
  ShapeType shape = ShapeType::DISC;
  float radius;
};

// without transform a ring is formed by a circle in the xy0 plane (not filled)
struct Ring : Shape
{
  ShapeType shape = ShapeType::RING;
  float radius;
};

// Suppose we have an untransformed cone with a certain radius and height
// We can truncate the cone by 'cutting' it open using two planes at z=a and z=b
// This results in 4 options:
//  - a and b both miss -> regular cone
//  - bottom is sliced off -> truncate
//  - top is sliced off -> truncate
//  - both are sliced off -> truncate
// For every case where we truncate we can either fill or not fill the resulting
// cut i.e. ring/disc
struct TruncatedCone : Shape
{
  ShapeType shape = ShapeType::TRUNCATED_CONE;
  float radius_bot, radius_top;
  float height;
  float bot, top; // if bot \in (0, height) a cut occurs
  bool fill_bot, fill_top;
};

// Think of a regular cone, but move the tip off-axis
// In this case, untransformed:
//   - tip = height*(cos(alpha), sin(alpha), 0) with alpha in [0, 90] degrees
struct ObliqueCone : Shape
{
  ShapeType shape = ShapeType::OBLIQUE_CONE;
  float radius;
  float height;
  float alpha; // in degrees
};
} // namespace Shapes

namespace DefaultFactories
{
class CircleLUT
{
  public:
  static float* LUTcosf;
  static float* LUTsinf;
  static int last_resolution;
  static void update(int resolution);

  private:
  CircleLUT()
  {
  }
  static void recompute(int resolution);
};

template <typename MeshType>
class SphereFactory
    : AbstractPrimitiveFactory<SphereFactory<MeshType>,
                               Primitives::Shapes::Sphere, MeshType>
{
  public:
  static MeshType createPrimitive(Primitives::Shapes::Sphere& shape,
                                  int resolution);
};

template <typename MeshType>
class RingFactory : AbstractPrimitiveFactory<RingFactory<MeshType>,
                                             Primitives::Shapes::Ring, MeshType>
{
  public:
  static MeshType createPrimitive(Primitives::Shapes::Ring& shape,
                                  int resolution);
};

template <typename MeshType>
class TruncatedConeFactory
    : AbstractPrimitiveFactory<TruncatedConeFactory<MeshType>,
                               Primitives::Shapes::TruncatedCone, MeshType>
{
  public:
  static MeshType createPrimitive(Primitives::Shapes::TruncatedCone& shape,
                                  int resolution);

  private:
  static MeshType last_ring;
  static int last_resolution;
};

} // namespace DefaultFactories
} // namespace Primitives

#endif