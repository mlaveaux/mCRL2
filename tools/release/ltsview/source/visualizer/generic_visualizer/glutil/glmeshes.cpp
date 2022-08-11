#include "glmeshes.h"
#include <fstream>
char space = ' ';
void GlUtil::MeshTypes::TriangleMesh::toObj(std::string filename){
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