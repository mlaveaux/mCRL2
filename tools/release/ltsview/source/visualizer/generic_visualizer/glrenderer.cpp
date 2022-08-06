#include "glrenderer.h"
#include "glscenegraph.h"


void Renderer::initialize(){
    if (!m_initialized){
        intializeShaders();
        m_initialized = true;
    }
}

void Renderer::forceInitialize(){
    intializeShaders();
    m_initialized = true;
}