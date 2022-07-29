#include "glrenderer.h"
#include "glscenegraph.h"
template<typename NodeData, typename ModelData>
void Renderer<NodeData, ModelData>::initialize(SceneGraph<NodeData, ModelData>* scenegraph){
    if (!m_initialized){
        intializeShaders();
        initializeData(scenegraph);
        m_initialized = true;
    }
}

template<typename NodeData, typename ModelData>
void Renderer<NodeData, ModelData>::forceInitialize(SceneGraph<NodeData, ModelData>* scenegraph){
    intializeShaders();
    initializeData(scenegraph);
    m_initialized = true;
}