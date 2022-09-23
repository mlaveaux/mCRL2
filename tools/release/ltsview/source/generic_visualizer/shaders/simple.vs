#version 330

uniform mat4 g_worldViewProjMatrix;
uniform float u_alpha = 1;

layout(location = 0) in vec3 v_vertex;
layout(location = 1) in vec3 v_normal;
layout(location = 2) in vec3 v_color;

out vec4 ex_color;

void main(void)
{
    gl_Position = g_worldViewProjMatrix * vec4(v_vertex, 1.0f);
    ex_color = vec4(v_color, u_alpha);
}