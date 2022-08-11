#pragma once
#ifndef MCRL2_OPENGL_MESHES_H
#define MCRL2_OPENGL_MESHES_H

#include <vector> 
#include <array>
#include <QVector3D>

namespace GlUtil
{
  
using Triangle = std::array<int, 4>;
using Quad = std::array<int, 4>;

namespace MeshTypes
{

struct Vertices{
  int n_vertices;
  std::vector<QVector3D> vertices;
};

struct TriangleMesh
{
  int n_vertices;
  std::vector<QVector3D> vertices;
  std::vector<QVector3D> vertex_normals;
  int n_triangles;
  std::vector<Triangle> triangles;
  std::vector<QVector3D> face_normals;

  void toObj(std::string filename);
};

struct TriangleMeshNoNormals
{
  int n_vertices;
  std::vector<QVector3D> vertices;
  int n_triangles;
  std::vector<Triangle> triangles;
};

struct QuadMesh
{
  int n_vertices;
  std::vector<QVector3D> vertices;
  std::vector<QVector3D> vertex_normals;
  int n_quads;
  std::vector<Quad> quads;
  std::vector<QVector3D> face_normals;
};

struct QuadMeshNoNormals
{
  int n_vertices;
  std::vector<QVector3D> vertices;
  int n_quads;
  std::vector<Quad> quads;
};
} // namespace MeshTypes
} // namespace GlUtil
#endif