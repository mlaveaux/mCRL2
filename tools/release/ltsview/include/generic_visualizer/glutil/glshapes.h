#ifndef MCRL2_OPENGL_GLSHAPES_H
#define MCRL2_OPENGL_GLSHAPES_H

namespace GlUtil
{
enum class ShapeType
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

struct Shape
{
  virtual ShapeType getShapeType() = 0;
};

namespace Shapes
{
struct Sphere : Shape
{
  ShapeType getShapeType() override
  {
    return ShapeType::SPHERE;
  }
  float radius = 1;
};

struct HemiSphere : Shape
{
  ShapeType getShapeType() override
  {
    return ShapeType::HEMISPHERE;
  }
  float radius;
};

// without transform a tube is formed by circles in the xy0 and xy1 planes
struct Tube : Shape
{
  ShapeType getShapeType() override
  {
    return ShapeType::TUBE;
  }
  float radius;
  float height;
};

// without transform a cone is formed by a cirle in the xy0 plane and a point at
// 001
struct Cone : Shape
{
  ShapeType getShapeType() override
  {
    return ShapeType::CONE;
  }
  float radius;
  float height;
};

// without transform a disc is formed by a circle in the xy0 plane (filled)
struct Disc : Shape
{
  ShapeType getShapeType() override
  {
    return ShapeType::DISC;
  }
  float radius;
};

// without transform a ring is formed by a circle in the xy0 plane (not filled)
struct Ring : Shape
{
  ShapeType getShapeType() override
  {
    return ShapeType::RING;
  }
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
  ShapeType getShapeType() override
  {
    return ShapeType::TRUNCATED_CONE;
  }
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
  ShapeType getShapeType() override
  {
    return ShapeType::OBLIQUE_CONE;
  }
  float radius;
  float height;
  float alpha; // in degrees
};
} // namespace Shapes
} // namespace GlUtil
#endif