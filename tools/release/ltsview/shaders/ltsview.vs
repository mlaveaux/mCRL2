#version 330
// per vertex
layout(location = 0) in vec3 vertex;
layout(location = 1) in vec3 normal;

// per shape
layout(location = 2) in vec4 color;
layout(location = 3) in mat4 model_mat;


uniform mat4 u_view;
uniform mat4 u_proj;
uniform float u_alpha;

out vec4 vs_color;

void main(void)
{
    vs_color = vec4(color.rgb, u_alpha);
    gl_Position = u_proj * u_view * model_mat * vec4(vertex, 1);
}