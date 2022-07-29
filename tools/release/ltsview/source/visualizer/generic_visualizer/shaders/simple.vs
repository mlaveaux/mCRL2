#version 330

uniform mat4 g_worldViewProjMatrix;

layout(location = 0) in vec3 vertex;

void main(void)
{
    gl_Position = g_worldViewProjMatrix * vec4(vertex, 1.0f);
}