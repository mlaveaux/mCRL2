#version 330

uniform vec4 g_color = vec4(1.0f);

out vec4 fragColor;

void main(void)
{
    fragColor = g_color;
}