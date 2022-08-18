#version 430


layout(std430, binding = 0) restrict readonly buffer b_vertices {
	vec4 vertices[];
};

layout(std430, binding = 1) restrict readonly buffer b_normals{
    vec4 normals[];
};

layout(std430, binding = 2) restrict readonly buffer b_colors {
    vec3 colors[];
};

layout(std430, binding = 3) restrict readonly buffer b_matrices {
    mat4 matrices[];
};

uniform mat4 u_view;
uniform mat4 u_proj;
uniform float u_alpha;

// struct VertexData{
//     int model_index;
//     int vertex_index;
// };
attribute ivec2 in_vertexData; 

out vec4 out_color;

void main(void)
{
    gl_Position =  u_proj * u_view * matrices[in_vertexData.x] * vertices[in_vertexData.y];
    out_color = vec4(colors[in_vertexData.x], 1);
}