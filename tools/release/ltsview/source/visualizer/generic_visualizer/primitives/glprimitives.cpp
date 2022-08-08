#include "glprimitives.h"
#include <cmath>
#include "mathutils.h"
using namespace MathUtils;


template < typename T, typename Shape, typename ModelData >
ModelData AbstractPrimitiveFactory<T, Shape, ModelData>::createPrimitive(Shape& shape, int resolution){
    return T::createPrimitive(shape, resolution);
}

using namespace Primitives;
using namespace Primitives::DefaultFactories;

// TODO: implement sphere
TriangleMesh SphereFactory<TriangleMesh>::createPrimitive(Shapes::Sphere& sphere, int resolution){
    return TriangleMesh();
}

// TODO: implement sphere
QuadMesh SphereFactory<QuadMesh>::createPrimitive(Shapes::Sphere& sphere, int resolution){
    return QuadMesh();
}

/// \brief Creates vertices for a ring. No normals/tris since ring does not define those
TriangleMesh RingFactory<TriangleMesh>::createPrimitive(Shapes::Ring& ring, int resolution){
    CircleLUT::update(resolution);
    int n_vertices = resolution;
    float* vertices = (float*) std::malloc(sizeof(float) * n_vertices * 3);
    int vertex_index = 0;
    float cosa, sina, x, y, z = 0;
    for (int i = 0; i < resolution; i++){
        cosa = CircleLUT::LUTcosf[i];
        sina = CircleLUT::LUTsinf[i];

        x = cosa * ring.radius;
        y = sina * ring.radius;
        vertices[vertex_index++] = x;
        vertices[vertex_index++] = y;
        vertices[vertex_index++] = z;
    }
    return {n_vertices, vertices, nullptr, 0, nullptr, nullptr};
}
float* CircleLUT::LUTcosf = nullptr;
float* CircleLUT::LUTsinf = nullptr;
int CircleLUT::last_resolution = -1;

void CircleLUT::update(int resolution){
    if (last_resolution != resolution) recompute(resolution);
}

void CircleLUT::recompute(int resolution){
    free(LUTcosf); free(LUTsinf); // clean arrays
    float* LUTcosf = (float*) std::malloc(sizeof(float) * resolution); // allocate space for arrays
    float* LUTsinf = (float*) std::malloc(sizeof(float) * resolution);
    float a = 0, da = (2.0*PI)/resolution;
    for (int i = 0; i < resolution; ++i) { LUTcosf[i] = cosf(a); LUTsinf[i] = sinf(a); a += da; }

}

int TruncatedConeFactory<TriangleMesh>::last_resolution = -1;
/// TODO: Look into possibly redefining into just the skewed cylinder
/// \brief Creates truncated cone out of triangles with bottom and/or top or neither closed
/// \param cone The cone to be created
/// \param resolution The number of vertices per circle. Resolution >= 3 required for normal behaviour
TriangleMesh TruncatedConeFactory<TriangleMesh>::last_ring; // we need to declare it, otherwise compiler doesn't know about it
TriangleMesh TruncatedConeFactory<TriangleMesh>::createPrimitive(Shapes::TruncatedCone& cone, int resolution){
    CircleLUT::update(resolution);

    int n_vertices = 2 * resolution;
    int n_triangles = n_vertices;
    if (cone.fill_bot) { n_triangles += resolution - 2; n_vertices += resolution; }
    if (cone.fill_top) { n_triangles += resolution - 2; n_vertices += resolution; }
    float* vertices       = (float*) std::malloc(sizeof(float) * n_vertices  * 3);
    float* vertex_normals = (float*) std::malloc(sizeof(float) * n_vertices  * 3);
    int*   triangles      = (int*)   std::malloc(sizeof(int)   * n_triangles * 3);
    float* face_normals   = (float*) std::malloc(sizeof(float) * n_triangles * 3);
    // create the vertices
    float cosa, sina, x, y, z = 0;
    int vertex_index = 0;
    int vertex_normal_index = 0;
    int triangle_index = 0;
    int face_normal_index = 0;
    float dr = cone.radius_bot - cone.radius_top;
    float inv_normal_magnitude = 1.0/sqrt(1+dr*dr);
    //bottom circle
    if (last_resolution != resolution){
        Shapes::Ring ring;
        ring.radius = 1;
        TruncatedConeFactory<TriangleMesh>::last_ring = RingFactory<TriangleMesh>::createPrimitive(ring, resolution);
        last_resolution = resolution;
    }
    float* bot_vertices = vertices;
    std::memcpy(vertices, last_ring.vertices, resolution*3*sizeof(float));
    std::transform(vertices, vertices+(resolution*3), vertices, [&cone](auto& c){ return c * cone.radius_bot; });
    
    float* top_vertices = vertices + (resolution*3);
    std::memcpy(top_vertices, last_ring.vertices, resolution*3*sizeof(float));
    std::transform(vertices, vertices+(resolution*3), vertices, [&cone](auto& c){ return c * cone.radius_top; });

    for (int i = 0; i < resolution; i++){
        cosa = CircleLUT::LUTcosf[i];
        sina = CircleLUT::LUTsinf[i];
        vertex_normals[vertex_normal_index++] = inv_normal_magnitude * cosa;
        vertex_normals[vertex_normal_index++] = inv_normal_magnitude * sina;
        vertex_normals[vertex_normal_index++] = inv_normal_magnitude * dr;
    }
    // for both cirlces the normals are identical
    std::memcpy(vertex_normals+(resolution*3), vertex_normals, resolution * 3 * sizeof(float));

    // triangulation of the sides
    int triangle_vertex_offset = resolution; // vertex number i of top circle is at location triangle_vertex_offset + i
    float cosa_next, sina_next;
    for (int i = 0; i < resolution; i++){
        int next = (i+1) % resolution;

        triangles[triangle_index++] = i;
        triangles[triangle_index++] = next;
        triangles[triangle_index++] = triangle_vertex_offset + i;

        triangles[triangle_index++] = next;
        triangles[triangle_index++] = triangle_vertex_offset + next;
        triangles[triangle_index++] = triangle_vertex_offset + i;

        // face normals
        cosa = CircleLUT::LUTcosf[i];
        sina = CircleLUT::LUTsinf[i];
        cosa_next = CircleLUT::LUTcosf[next];
        sina_next = CircleLUT::LUTsinf[next];

        face_normals[face_normal_index++] = inv_normal_magnitude * ONE_THIRD * (2*cosa + cosa_next);
        face_normals[face_normal_index++] = inv_normal_magnitude * ONE_THIRD * (2*sina + sina_next);
        face_normals[face_normal_index++] = inv_normal_magnitude * dr;

        face_normals[face_normal_index++] = inv_normal_magnitude * ONE_THIRD * (cosa + 2*cosa_next);
        face_normals[face_normal_index++] = inv_normal_magnitude * ONE_THIRD * (sina + 2*sina_next);
        face_normals[face_normal_index++] = inv_normal_magnitude * dr;
    }

    if (cone.fill_bot || cone.fill_top){
        float* top_disk = vertices + (2*resolution*3);
        std::memcpy(top_disk, top_vertices, resolution*3*sizeof(float));
        float* bot_disk = top_disk + (cone.fill_top ? resolution*3 : 0);
        std::memcpy(bot_disk, bot_vertices, resolution*3*sizeof(float));
        float* top_vertex_normals = vertex_normals + (2*resolution*3);
        float* bot_vertex_normals = top_vertex_normals + (cone.fill_top ? resolution*3 : 0);
        
        for (int i = 0; i < resolution; i++){
            top_vertex_normals[i*3+2] = 1;
            bot_vertex_normals[i*3+2] = -1;
        }
        if (cone.fill_bot){
            // triangulation of the circle iff cone.fill_bot by connecting:
            // 0->2->1, 0->3->2 etc. ~ 0->i->i-1
            // Reason: OpenGL by default uses CCW winding
            for (int i = 2; i < resolution; i++){
                triangles[triangle_index++] = 0;
                triangles[triangle_index++] = i;
                triangles[triangle_index++] = i-1;

                // The normal for a given triangle on the bottom is always straight down
                face_normals[face_normal_index++] = 0;
                face_normals[face_normal_index++] = 0;
                face_normals[face_normal_index++] = -1;
            }
        }

        if (cone.fill_top){
            // top circle triangulation iff cone.fill_top
            // triangulation now goes 0-1-2 with 0-1-2 offset by number of vertices
            for (int i = 2; i < n_triangles && cone.fill_top; i++){
                triangles[triangle_index++] = triangle_vertex_offset+0;
                triangles[triangle_index++] = triangle_vertex_offset+i-1;
                triangles[triangle_index++] = triangle_vertex_offset+i;

                // The normal for a given triangle on the top is always straight up
                face_normals[face_normal_index++] = 0;
                face_normals[face_normal_index++] = 0;
                face_normals[face_normal_index++] = 1;
            }
        }
    }
    return {n_vertices, vertices, vertex_normals, n_triangles, triangles, face_normals};
}

// TODO: implement truncated cone
// QuadMesh DefaultFactories::TruncatedConeFactory<QuadMesh>::createPrimitive(Shapes::TruncatedCone& cone, int resolution){
//     return QuadMesh();
// }