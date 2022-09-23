#include "glprimitivefactories.h"
#include "mathutils.h"

using namespace GlUtil;

template <typename T, typename Shape, typename MeshType>
MeshType
AbstractPrimitiveFactory<T, Shape, MeshType>::createPrimitive(Shape* shape,
                                                              int resolution)
{
  return T::createPrimitive(shape, resolution);
}

using namespace DefaultFactories;

/// \brief Creates vertices for a ring. No normals/tris since ring does not
/// define those
template<>
Meshes::MeshTypes::Vertices
RingFactory<Meshes::MeshTypes::Vertices>::createPrimitive(Shapes::Ring* ring,
                                                  int resolution)
{
  LUTs::CircleLUT::update(resolution);
  int n_vertices = resolution;
  std::vector<QVector3D> vertices(n_vertices);
  int vertex_index = 0;
  qreal cosa, sina, x, y, z = 0;
  for (int i = 0; i < resolution; i++)
  {
    cosa = LUTs::CircleLUT::LUTcosf[i];
    sina = LUTs::CircleLUT::LUTsinf[i];

    x = cosa * ring->radius;
    y = sina * ring->radius;
    vertices[vertex_index++] = {(float)x, (float)y, (float)z};
  }
  return {n_vertices, vertices};
}

uint encode(int a, int b)
{
  return ((uint)std::min(a, b) << 16) + (uint)std::max(a, b);
}
/// \brief Creates a trianglemesh of an icosahedron subdivided \param resolution
/// times. \param resolution Number of subdivisions. Assumed to be <= 6 (~41k
/// verts ~82k tris)
// Credit for vertex/triangle indices:
// http://blog.andreaskahler.com/2009/06/creating-icosphere-mesh-in-code.html
template <>
Meshes::MeshTypes::TriangleMesh
SphereFactory<Meshes::MeshTypes::TriangleMesh>::createPrimitive(Shapes::Sphere* sphere,
                                                        int resolution)
{
  int n_vertices = 12;
  int n_triangles = 20;

  qreal sqrt5 = sqrt(5);
  float s = sqrtf((5 - sqrt5) * 0.1);
  float t = sqrtf((5 + sqrt5) * 0.1);

  std::vector<QVector3D> vertices = {
      {-s, t, 0},  {s, t, 0},  {-s, -t, 0}, {s, -t, 0}, {0, -s, t},  {0, s, t},
      {0, -s, -t}, {0, s, -t}, {t, 0, -s},  {t, 0, s},  {-t, 0, -s}, {-t, 0, s},
  };

  // all counter clockwise winding
  std::vector<GlUtil::Triangle> triangles = {
      // 5 faces around point 0
      {0, 11, 5},
      {0, 5, 1},
      {0, 1, 7},
      {0, 7, 10},
      {0, 10, 11},

      // 5 adjacent faces
      {1, 5, 9},
      {5, 11, 4},
      {11, 10, 2},
      {10, 7, 6},
      {7, 1, 8},

      // 5 faces around point 3
      {3, 9, 4},
      {3, 4, 2},
      {3, 2, 6},
      {3, 6, 8},
      {3, 8, 9},

      // 5 adjacent faces
      {4, 9, 5},
      {2, 4, 11},
      {6, 2, 10},
      {8, 6, 7},
      {9, 8, 1},
  };

  // now we have to subdivide.
  // every subdivision creates 1 points per edge;
  //  every face has 3 edges, but every edge is shared. and turns one face into
  //  4 faces.
  // In short:
  //   n_vertices_new = n_vertices + 3*(n_triangles)/2
  //   n_triangles_new = 4*n_triangles
  while (resolution-- > 0)
  {
    std::map<uint, int> midpointmap = std::map<uint, int>();
    // update counts
    int a, b, c;
    int ab, ac, bc;
    int triangle_index = 0;
    int vertex_index = n_vertices;
    n_vertices += 3 * (n_triangles >> 1);
    n_triangles *= 4;
    vertices.resize(n_vertices);
    std::vector<GlUtil::Triangle> new_triangles =
        std::vector<GlUtil::Triangle>(n_triangles);
    // Loop over triangles and subdivide
    for (int i = 0; i < n_triangles >> 2; i++)
    {
      a = triangles[i][0];
      b = triangles[i][1];
      c = triangles[i][2];

      // check if halfway points exist
      auto ab_it = midpointmap.find(encode(a, b));
      if (ab_it != midpointmap.end())
      {
        ab = ab_it->second;
      }
      else
      {
        ab = vertex_index;
        vertices[ab] = (vertices[a] + vertices[b]).normalized();
        midpointmap.insert({encode(a, b), ab});
        vertex_index++;
      }

      auto ac_it = midpointmap.find(encode(a, c));
      if (ac_it != midpointmap.end())
      {
        ac = ac_it->second;
      }
      else
      {
        ac = vertex_index;
        vertices[ac] = (vertices[a] + vertices[c]).normalized();
        midpointmap.insert({encode(a, c), ac});
        vertex_index++;
      }

      auto bc_it = midpointmap.find(encode(b, c));
      if (bc_it != midpointmap.end())
      {
        bc = bc_it->second;
      }
      else
      {
        bc = vertex_index;
        vertices[bc] = (vertices[b] + vertices[c]).normalized();
        midpointmap.insert({encode(b, c), bc});
        vertex_index++;
      }

      // since we know CCW winding, we know which triangles will be created
      new_triangles[triangle_index++] = {a, ab, ac};
      new_triangles[triangle_index++] = {b, bc, ab};
      new_triangles[triangle_index++] = {c, ac, bc};
      new_triangles[triangle_index++] = {ab, bc, ac};
    }
    triangles = new_triangles;
  }

  /// TODO: Point vertex_normals to vertices array.
  // for a unit sphere, any vertex points the same way as the normal
  std::vector<QVector3D> vertex_normals(vertices);

  // Face normals are the average of the associated 3 vertex normals
  std::vector<QVector3D> face_normals(triangles.size());
  int normal_index = 0;
  int a, b, c;
  for (int i = 0; i < n_triangles; i++)
  {
    a = triangles[i][0];
    b = triangles[i][1];
    c = triangles[i][2];
    face_normals[normal_index++] =
        (float)MathUtils::ONE_THIRD *
        (vertex_normals[a] + vertex_normals[b] + vertex_normals[c]);
  }

  // finally we need to multiply all vertices by the radius of the sphere
  qreal r = sphere->radius;
  for (auto it = vertices.begin(); it != vertices.end(); ++it)
  {
    *it *= r;
  }
  return {n_vertices,  vertices,  vertex_normals,
          n_triangles, triangles, face_normals};
}

/// TODO: implement quad sphere
template<>
Meshes::MeshTypes::QuadMesh
SphereFactory<Meshes::MeshTypes::QuadMesh>::createPrimitive(Shapes::Sphere* sphere,
                                                    int resolution)
{
  return Meshes::MeshTypes::QuadMesh();
}


template<>
Meshes::MeshTypes::Vertices TruncatedConeFactory<
    Meshes::MeshTypes::TriangleMesh>::last_ring = Meshes::MeshTypes::Vertices(); // we need to declare it, otherwise
                                         // compiler doesn't know about it
template<>
int TruncatedConeFactory<Meshes::MeshTypes::TriangleMesh>::last_resolution = -1;
/// TODO: Look into possibly redefining into just the skewed cylinder without
/// the top and bottom lids \brief Creates truncated cone out of triangles with
/// bottom and/or top or neither closed \param cone The cone to be created
/// \param resolution The number of vertices per circle. Resolution >= 3
/// required for normal behaviour
template<>
Meshes::MeshTypes::TriangleMesh
TruncatedConeFactory<Meshes::MeshTypes::TriangleMesh>::createPrimitive(
    Shapes::TruncatedCone* cone, int resolution)
{
  LUTs::CircleLUT::update(resolution);

  int n_vertices = 2 * resolution;
  int n_triangles = n_vertices;
  if (cone->fill_bot)
  {
    n_triangles += resolution - 2;
    n_vertices += resolution;
  }
  if (cone->fill_top)
  {
    n_triangles += resolution - 2;
    n_vertices += resolution;
  }
  std::vector<QVector3D> vertices(n_vertices);
  std::vector<QVector3D> vertex_normals(n_vertices);
  std::vector<GlUtil::Triangle> triangles(n_triangles);
  std::vector<QVector3D> face_normals(n_triangles);
  // create the vertices
  qreal cosa, sina, z = 0;
  int vertex_index = 0;
  int vertex_normal_index = 0;
  int triangle_index = 0;
  int face_normal_index = 0;
  qreal dr = cone->radius_bot - cone->radius_top;
  qreal inv_normal_magnitude = 1.0 / sqrt(1 + dr * dr);
  // bottom circle
  if (last_resolution != resolution)
  {
    Shapes::Ring* ring = new Shapes::Ring();
    ring->radius = 1;
    TruncatedConeFactory<Meshes::MeshTypes::TriangleMesh>::last_ring =
        RingFactory<Meshes::MeshTypes::Vertices>::createPrimitive(ring, resolution);
    last_resolution = resolution;
  }
  std::copy(last_ring.vertices.begin(), last_ring.vertices.end(),
            vertices.begin());
  std::transform(vertices.begin(), vertices.begin() + resolution,
                 vertices.begin(),
                 [&cone](auto& c) { return c * cone->radius_bot; });

  std::copy(last_ring.vertices.begin(), last_ring.vertices.end(),
            vertices.begin() + resolution);
  QVector3D top(0, 0, 1);
  std::transform(vertices.begin() + resolution,
                 vertices.begin() + 2 * resolution,
                 vertices.begin() + resolution,
                 [&cone](auto& c) { return c * cone->radius_top; });
  for (int i = resolution; i < 2*resolution; i++){
    vertices[i][2] += cone->height;
  }

  for (int i = 0; i < resolution; i++)
  {
    cosa = LUTs::CircleLUT::LUTcosf[i];
    sina = LUTs::CircleLUT::LUTsinf[i];
    vertex_normals[vertex_normal_index++] =
        inv_normal_magnitude * QVector3D({(float)cosa, (float)sina, (float)dr});
  }
  // for both cirlces the normals are identical
  std::copy(vertex_normals.begin(), vertex_normals.begin() + resolution,
            vertex_normals.begin() + resolution);

  // triangulation of the sides
  int triangle_vertex_offset =
      resolution; // vertex number i of top circle is at location
                  // triangle_vertex_offset + i
  qreal cosa_next, sina_next;
  for (int i = 0; i < resolution; i++)
  {
    int next = (i + 1) % resolution;

    triangles[triangle_index++] = {i, next, triangle_vertex_offset + i};

    triangles[triangle_index++] = {next, triangle_vertex_offset + next,
                                   triangle_vertex_offset + i};

    // face normals
    cosa = LUTs::CircleLUT::LUTcosf[i];
    sina = LUTs::CircleLUT::LUTsinf[i];
    cosa_next = LUTs::CircleLUT::LUTcosf[next];
    sina_next = LUTs::CircleLUT::LUTsinf[next];

    face_normals[face_normal_index++] = {
        (float)(inv_normal_magnitude * MathUtils::ONE_THIRD * (2 * cosa + cosa_next)),
        (float)(inv_normal_magnitude * MathUtils::ONE_THIRD * (2 * sina + sina_next)),
        (float)(inv_normal_magnitude * dr)};

    face_normals[face_normal_index++] = {
        (float)(inv_normal_magnitude * MathUtils::ONE_THIRD * (cosa + 2 * cosa_next)),
        (float)(inv_normal_magnitude * MathUtils::ONE_THIRD * (sina + 2 * sina_next)),
        (float)(inv_normal_magnitude * dr)};
  }

  if (cone->fill_bot || cone->fill_top)
  {
    auto bot_vertices = vertices.begin();
    auto bot_disk = vertices.begin() + (2 * resolution);
    std::copy(bot_vertices, bot_vertices + resolution, bot_disk);
    auto top_vertices = vertices.begin() + resolution;
    auto top_disk = bot_disk + (cone->fill_bot ? resolution : 0);
    if (cone->fill_top){
      std::copy(top_vertices, top_vertices + resolution, top_disk);
    }

    auto bot_vertex_normals = vertex_normals.begin() + (2 * resolution);
    auto top_vertex_normals =
        bot_vertex_normals + (cone->fill_bot ? resolution : 0);

    for (int i = 0; i < resolution; i++)
    {
      if (cone->fill_bot)
        bot_vertex_normals[i][2] = -1;
      if (cone->fill_top)
        top_vertex_normals[i][2] = 1;
    }
    int vertex_offset = 2*resolution;

    if (cone->fill_bot)
    {
      // triangulation of the circle iff cone.fill_bot by connecting:
      // 0->2->1, 0->3->2 etc. ~ 0->i->i-1
      // Reason: OpenGL by default uses CCW winding
      for (int i = 2; i < resolution; i++)
      {
        triangles[triangle_index++] = {vertex_offset+0, vertex_offset+i, vertex_offset+i - 1};

        // The normal for a given triangle on the bottom is always straight down
        face_normals[face_normal_index++] = {0, 0, -1};
      }
      vertex_offset += resolution;
    }
    
    if (cone->fill_top)
    {
      // top circle triangulation iff cone.fill_top
      // triangulation now goes 0-1-2 with 0-1-2 offset by number of vertices
      for (int i = 2; i < resolution; i++)
      {
        triangles[triangle_index++] = {vertex_offset + 0,
                                       vertex_offset + i - 1,
                                       vertex_offset + i};

        // The normal for a given triangle on the top is always straight up
        face_normals[face_normal_index++] = {0, 0, 1};
      }
    }
  }
  return {n_vertices,  vertices,  vertex_normals,
          n_triangles, triangles, face_normals};
}