#version 430 core


layout(std140) buffer b_vertices {
	vec4 vertices[];
};

layout(std140) buffer b_normals{
    vec4 normals[];
};

layout(std140) buffer b_colors {
    vec4 colors[];
};

layout(std140) buffer b_matrices {
    mat4 matrices[];
};

uniform mat4 u_view;
uniform mat4 u_proj;
uniform float u_alpha;
// struct VertexData{
//     int model_index;
//     int vertex_index;
// };
in ivec2 in_vertexData; 

out vec4 vs_color;

void main(void)
{
    // gl_Position =  u_proj * u_view * matrices[in_vertexData.x] * vertices[in_vertexData.y];
    int x = in_vertexData.x + in_vertexData.y;
    vs_color = u_proj * u_view * vec4(x, x, x, u_alpha);
    gl_Position = vec4(x);
}