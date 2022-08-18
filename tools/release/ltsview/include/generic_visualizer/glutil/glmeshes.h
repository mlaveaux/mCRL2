#pragma once
#ifndef MCRL2_OPENGL_MESHES_H
#define MCRL2_OPENGL_MESHES_H
#define INF 10000000000000

#include <vector>
#include <array>
#include <QVector3D>
#include <fstream>
#include <sstream>
#include <cmath>

namespace GlUtil
{

using Triangle = std::array<int, 3>;
using Quad = std::array<int, 4>;
namespace Meshes
{

struct AABB{
  QVector3D min, max;
  QVector3D getBounds(){ return max-min; }
};

namespace MeshTypes
{

struct Vertices
{
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


template < typename MeshType >
AABB toAABB(MeshType mesh){
  QVector3D min = QVector3D(INF, INF, INF);
  QVector3D max = QVector3D(-INF, -INF, -INF);
  for (auto& vert : mesh.vertices){
    min.setX(std::min(min.x(), vert.x()));
    min.setY(std::min(min.y(), vert.y()));
    min.setZ(std::min(min.z(), vert.z()));
    max.setX(std::max(max.x(), vert.x()));
    max.setY(std::max(max.y(), vert.y()));
    max.setZ(std::max(max.z(), vert.z()));
  }
  return {min, max};
}

namespace Export{
  enum class FileType{
    PLY,
    OBJ
  };
  template <typename MeshType, FileType ftype>
  void exportMesh(std::string& filename, MeshType& mesh);
  
  template <>
  void exportMesh<MeshTypes::TriangleMesh, FileType::PLY>(std::string& filename, MeshTypes::TriangleMesh& mesh);
}
} // namespace Meshes
} // namespace GlUtil
#endif