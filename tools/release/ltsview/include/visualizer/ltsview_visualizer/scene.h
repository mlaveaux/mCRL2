#ifndef SCENE_H
#define SCENE_H

#include "glscene.h"
#include <vector>
#include "cluster.h"

namespace Scene{
    struct NodeData
    {
        int model_id;
        int model_matrix_id;
    };

    struct Model
    {
        int n_vertices;
        float* vertices;
        int n_triangles;
        int* triangles;
    };

    struct ModelData{
        std::vector<Model> models;
        std::vector<QMatrix4x4> matrices;
    };
    
    struct SceneGraphFunctor{
        
    };
}



class Scene : GLScene<Scene::NodeData, Scene::ModelData>{
  public:
  void generateSceneGraph(Cluster* root);
};

#endif