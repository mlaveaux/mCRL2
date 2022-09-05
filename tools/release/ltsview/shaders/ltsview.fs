#version 430 core

in vec4 vs_color;

void main(void)
{
  gl_FragColor = vs_color;
}