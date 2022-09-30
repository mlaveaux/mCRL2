#version 430 compatibility

in vec4 vs_color;
out vec4 out_color;

void main(void)
{
  out_color = vs_color;
}