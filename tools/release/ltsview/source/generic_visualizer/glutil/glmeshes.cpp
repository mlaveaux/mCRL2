#include "glmeshes.h"
#include <fstream>
#include <cassert>
char space = ' ';
/// TODO: Move to GlUtil::Meshes::Exporters
void GlUtil::Meshes::MeshTypes::TriangleMesh::toObj(std::string filename){
    std::ofstream file(filename);
    if (!file) return;
    
    for (auto& vert : vertices){
        file << "v " << vert[0] << space << vert[1] << space << vert[2] << std::endl;
    }

    for (auto& norm : vertex_normals){
        file << "vn " << norm[0] << space << norm[1] << space << norm[2] << std::endl;
    }

    for (auto& tri : triangles){
        file << "f " << tri[0]+1 << "//" << tri[0]+1 << space << tri[1]+1 << "//" << tri[1]+1 << space << tri[2]+1 << "//" << tri[2]+1 << std::endl;
    }

    file.close();
}

template <>
void GlUtil::Meshes::Export::exportMesh<GlUtil::Meshes::MeshTypes::TriangleMesh, GlUtil::Meshes::Export::FileType::PLY>(std::string& filename, GlUtil::Meshes::MeshTypes::TriangleMesh& mesh){
    std::ofstream file(filename);
    std::stringstream string_file;
    assert(file);

    if (!file)
    {
      return ;
    }

    string_file << "ply" << std::endl
         << "format ascii 1.0" << std::endl
         << "element vertex " << mesh.n_vertices << std::endl;
    string_file << "property float x" << std::endl;
    string_file << "property float y" << std::endl;
    string_file << "property float z" << std::endl;
    string_file << "property float nx" << std::endl;
    string_file << "property float ny" << std::endl;
    string_file << "property float nz" << std::endl;
    string_file << "element face " << mesh.n_triangles << std::endl;
    string_file << "property list uchar uint vertex_indices" << std::endl;
    string_file << "end_header" << std::endl;
    
    // dump ply file
    char space = ' ';
    for (int i = 0; i < mesh.n_vertices; i++){
        auto& vert = mesh.vertices[i];
        auto& norm = mesh.vertex_normals[i];
        string_file << vert.x() << space << vert.y() << space << vert.z() << space
             << norm.x() << space << norm.y() << space << norm.z()
             << std::endl;
    }
    for (auto& tri : mesh.triangles){
      string_file <<  "3 " <<  tri[0] << space << tri[1] << space << tri[2] << std::endl;
    }
    file << string_file.rdbuf();
    file.close();
  }